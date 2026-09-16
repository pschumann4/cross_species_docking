"""
mds_structure_prep.py
=============================
Prepares the test and reference structure for molecular dynamics simulation and flexibility analysis. 
It performs the following steps:
1. Runs PDBFixer to clean up the structure, add missing atoms, and ensure it's suitable for MD.
2. Identifies binding pocket residues based on proximity to the ligand (in reference structure).
3. Runs a fixed-length (MDS_TIME_PS) MD simulation, seeded for reproducibility, to relax the
   predicted model and sample within-basin conformational fluctuations.
4. Analyzes the trajectory to calculate RMSF (Root Mean Square Fluctuation) for each residue,
   identifying which residues are flexible based on a defined threshold.
5. Extracts an ensemble of N_ENSEMBLE cluster medoids (by pocket RMSD) from the post-burn-in
   frames, capturing diverse relaxed conformations to dock against.

This MD step is for RELAXATION + FLEXIBILITY PROFILING + within-basin ensemble sampling. It does
NOT detect thermodynamic equilibrium and does NOT reach distinct conformational macrostates on this
timescale.

Outputs
-------
- Fixed PDB files for test and reference structures (protein only for MD, with ligand for reference)
- A relaxation QC plot (superposed Cα RMSD vs first frame) per structure
- RMSF plots showing residue flexibility with flexible residues highlighted
- A summary file listing flexible residues for each structure
- N_ENSEMBLE medoid ensemble structures per receptor (with intra-ensemble diversity QC)
"""

import os
import sys
import shutil
import time
import warnings
import numpy as np

# MDAnalysis.analysis.align pulls in the deprecated Bio.Application module, which
# emits a BiopythonDeprecationWarning at import time. It is a benign third-party
# warning (we do not use Bio.Application); silence it before importing MDAnalysis
# so it does not clutter the pipeline output.
warnings.filterwarnings("ignore", message=".*Bio.Application modules.*")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import euclidean3d, centroid, load_config, save_config, CONFIG_FILENAME, prompt_yes_no
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis import align
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# ---------------------------------------------------------------------------
# MD protocol constants (fixed for reproducibility; no longer user-tunable)
# ---------------------------------------------------------------------------
MDS_TIME_PS = 200          # fixed simulation length (relaxation + sampling)
INTEGRATOR_SEED = 42       # Langevin seed → reproducible trajectories per machine
BURN_IN_FRACTION = 0.25    # leading fraction of frames discarded as relaxation transient
N_ENSEMBLE = 5             # number of medoid structures extracted per receptor
LOW_DIVERSITY_RMSD = 0.5   # Å; below this mean pairwise medoid RMSD the pocket is
                           # effectively rigid (flagged as a QC warning)


# ---------------------------------------------------------------------------
# Ensemble-selection helpers (pure; unit-tested in tests/test_utils.py)
# ---------------------------------------------------------------------------

def apply_burn_in(n_frames, fraction=BURN_IN_FRACTION):
    """
    Return the index of the first frame to keep after discarding a leading
    burn-in fraction (the relaxation transient where the structure is still
    shedding the predicted-model geometry).

    Guarantees at least two frames remain so a pairwise distance matrix can be
    built; for very short trajectories the burn-in is reduced accordingly.
    """
    if n_frames < 2:
        return 0
    start = int(n_frames * fraction)
    return min(start, n_frames - 2)


def select_medoids(dist_matrix, k):
    """
    Cluster frames by a precomputed pairwise distance matrix and return one
    medoid index per cluster.

    Hierarchical (average-linkage) clustering partitions the frames into at most
    k clusters; each cluster's medoid is the member minimizing the summed
    distance to the other members — i.e. a real, representative frame. Returns
    sorted local indices into dist_matrix.

    When there are k or fewer frames, every frame is returned (nothing to merge).
    """
    n = len(dist_matrix)
    if n == 0:
        return []
    if n <= k:
        return list(range(n))

    condensed = squareform(np.asarray(dist_matrix), checks=False)
    linkage_matrix = linkage(condensed, method="average")
    labels = fcluster(linkage_matrix, t=k, criterion="maxclust")

    medoids = []
    for label in np.unique(labels):
        members = np.where(labels == label)[0]
        # medoid = member with the smallest total distance to the rest of its cluster
        submatrix = np.asarray(dist_matrix)[np.ix_(members, members)]
        medoid_local = members[np.argmin(submatrix.sum(axis=1))]
        medoids.append(int(medoid_local))
    return sorted(medoids)


def intra_ensemble_diversity(dist_matrix, medoid_indices):
    """
    Mean pairwise distance among the selected medoids — a QC measure of how much
    conformational diversity the ensemble actually captures. Returns 0.0 when
    fewer than two medoids are present.
    """
    if len(medoid_indices) < 2:
        return 0.0
    dm = np.asarray(dist_matrix)
    pairs = [
        dm[medoid_indices[i], medoid_indices[j]]
        for i in range(len(medoid_indices))
        for j in range(i + 1, len(medoid_indices))
    ]
    return float(np.mean(pairs))


def run_pdbfixer(pdb_name, output_dir, keep_ligand=False, ligand_name=None):
    """
    Run PDBFixer on a PDB file to fix common issues
    
    Parameters:
    -----------
    pdb_name : str
        Path to the input PDB file
    output_dir : str
        Directory to save output files
    keep_ligand : bool
        Whether to keep heteroatoms (ligands) in the structure
    ligand_name : str, optional
        Specific ligand residue name to keep (e.g., 'BNF')
    
    Returns:
    --------
    tuple: (fixed_filename, original_resids)
        - fixed_filename: Path to the fixed PDB file (protein only for MD, or with ligand for reference)
        - original_resids: Array of original residue IDs to preserve numbering
    """
    os.makedirs(output_dir, exist_ok=True)
    
    print("\nPreparing", os.path.basename(pdb_name), "using PDBFixer...")
    
    # Store original residue IDs before any processing
    original_u = mda.Universe(pdb_name)
    original_resids = original_u.select_atoms("protein").residues.resids.copy()
    del original_u
    
    # Fix the PDB file
    fixer = PDBFixer(pdb_name)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    
    # Always remove heterogens for the fixed file (needed for MD)
    fixer.removeHeterogens(keepWater=False)
    
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    # Output fixed PDB (protein only)
    fixed_filename = os.path.join(output_dir, os.path.basename(pdb_name).replace(".pdb", "_fixed.pdb"))
    PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_filename, 'w'))
    
    # If this is a reference structure with ligand, append ligand to create final version
    if keep_ligand and ligand_name:
        print(f"  Preserving ligand {ligand_name} from original structure...")
        
        # Read the original PDB to get ligand lines
        with open(pdb_name, 'r') as f:
            original_lines = f.readlines()
        
        # Extract ligand lines
        ligand_lines = [line for line in original_lines 
                       if line.startswith('HETATM') and line[17:20].strip() == ligand_name]
        
        if ligand_lines:
            # Append ligand to fixed file (before END record)
            with open(fixed_filename, 'r') as f:
                fixed_lines = f.readlines()
            
            # Find where to insert (before TER or END)
            insert_index = len(fixed_lines)
            for i, line in enumerate(fixed_lines):
                if line.startswith('TER') or line.startswith('END'):
                    insert_index = i
                    break
            
            # Insert ligand lines
            new_lines = fixed_lines[:insert_index] + ligand_lines + fixed_lines[insert_index:]
            
            with open(fixed_filename, 'w') as f:
                f.writelines(new_lines)
            
            print(f"  Ligand {ligand_name} preserved ({len(ligand_lines)} atoms)")
        else:
            print(f"  Warning: Ligand {ligand_name} not found in original file")

    return fixed_filename, original_resids


# Solvent residue names that are never the ligand of interest; excluded from the
# candidate ligand list unless they are the only HETATM records present.
WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SPC"}


def list_hetatm_ligands(pdb_file):
    """
    Return candidate ligand residue names found in a PDB's HETATM records.

    Returns a list of (resname, atom_count) sorted by total atom count
    descending — the true ligand is usually the largest HETATM group, while
    ions and waters are 1-3 atoms. Water-like residues are dropped unless they
    are the only HETATM records present.
    """
    counts = {}
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("HETATM"):
                    resname = line[17:20].strip()
                    if resname:
                        counts[resname] = counts.get(resname, 0) + 1
    except OSError as e:
        print(f"  Warning: could not read {pdb_file}: {e}")
        return []

    filtered = {k: v for k, v in counts.items() if k.upper() not in WATER_RESNAMES}
    if not filtered:                      # only water present — show it anyway
        filtered = counts
    return sorted(filtered.items(), key=lambda kv: kv[1], reverse=True)


def get_binding_pocket_residues(pdb_file, ligand_name):
    """
    Identify binding pocket residues for RMSF analysis based on proximity to ligand
    
    Parameters:
    -----------
    pdb_file : str
        Path to PDB file
    ligand_name : str
        3-letter ligand code
    
    Returns:
    --------
    list: Residue IDs (as strings) of binding pocket residues
    """
    # Read PDB file
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    
    # Extract ligand coordinates
    ligand_coords = []
    for line in lines:
        if line.startswith("HETATM") and line[17:20].strip() == ligand_name:
            ligand_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    
    if not ligand_coords:
        print(f"Warning: Ligand {ligand_name} not found in {os.path.basename(pdb_file)}")
        return []
    
    # Calculate ligand centroid
    ligand_centroid = centroid(ligand_coords)
    
    # Extract residue coordinates
    res_coords = {}
    for line in lines:
        if line.startswith("ATOM"):
            resnr = line[22:26].strip()
            coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if resnr not in res_coords:
                res_coords[resnr] = [coords]
            else:
                res_coords[resnr].append(coords)
    
    # Calculate residue centroids
    res_centroids = {resnr: centroid(coords) for resnr, coords in res_coords.items()}
    
    # Define binding site distance
    BS_DIST = 7.5
    
    # Calculate ligand radius (max distance from centroid to any atom)
    max_dist = max(euclidean3d(ligand_centroid, coord) for coord in ligand_coords)
    cutoff = BS_DIST + max_dist
    
    # Find residues within cutoff based on centroid distance
    bindingsite_resnr = [resnr for resnr, res_cent in res_centroids.items() 
                         if euclidean3d(ligand_centroid, res_cent) < cutoff]
    
    # Also check minimum atom-atom distance < 7.5 Å…
    min_dist_residues = []
    for resnr, coords in res_coords.items():
        min_d = min(euclidean3d(lig_coord, res_coord) 
                   for lig_coord in ligand_coords 
                   for res_coord in coords)
        if min_d < 7.5:
            min_dist_residues.append(resnr)
    
    # Take intersection
    bindingsite_resnr = [r for r in bindingsite_resnr if r in min_dist_residues]
    
    return sorted(bindingsite_resnr)


def save_flexible_residues_summary(output_dir, results):
    """
    Create a consolidated flexible residues file
    
    Expected format:
    structure_name: [resid1,resid2,resid3]
    
    Parameters:
    -----------
    output_dir : str
        Output directory
    results : list
        List of result dictionaries from processing
    """
    os.makedirs(output_dir, exist_ok=True)
    flex_summary_path = os.path.join(output_dir, 'flex_residues.txt')
    
    # Collect structures with flexibility data
    flex_structures = []
    for result in results:
        if result.get('rmsf_data') and not result.get('is_reference', False):
            rmsf_data = result['rmsf_data']
            structure_name = result['structure_name']
            flexible_residues = rmsf_data['flexible_residues']
            
            flex_structures.append({
                'name': structure_name,
                'residues': flexible_residues,
                'count': len(flexible_residues)
            })
    
    if not flex_structures:
        print("\nNo flexible residue data to save.")
        return None
    
    # Calculate statistics using explicit builtin sum
    total_structures = len(flex_structures)
    structures_with_flex = len([s for s in flex_structures if s['count'] > 0])
    
    # Write summary file
    with open(flex_summary_path, 'w') as f:
        f.write("FLEXIBLE RESIDUES SUMMARY\n")
        f.write("="*60 + "\n\n")
        
        for struct in sorted(flex_structures, key=lambda x: x['name']):
            # Convert residue IDs to comma-separated list
            if len(struct['residues']) > 0:
                resid_list = ','.join(str(r) for r in sorted(struct['residues']))
                f.write(f"{struct['name']}: [{resid_list}]\n")
            else:
                f.write(f"{struct['name']}: []\n")
        
        f.write("\n" + "="*60 + "\n")
        f.write(f"Total structures analyzed: {total_structures}\n")
        f.write(f"Structures with flexible residues: {structures_with_flex}\n")
    
    print(f"\n  Flexible residues summary saved to: {flex_summary_path}")
    
    return flex_summary_path


def calculate_rmsf(trajectory_path, binding_pocket_resids=None, dt=5.0, original_resids=None):
    """
    Calculate RMSF (Root Mean Square Fluctuation) for protein residues
    
    Parameters:
    -----------
    trajectory_path : str
        Path to trajectory PDB file
    binding_pocket_resids : list, optional
        List of residue IDs (as strings) to focus on. If None, calculates for all residues
    dt : float
        Time between frames in ps
    original_resids : np.array, optional
        Original residue numbering to map back to
    
    Returns:
    --------
    dict: Contains 'resids', 'rmsf_values', 'flexible_residues' (RMSF > threshold)
    """
    u = None
    try:
        u = mda.Universe(trajectory_path, dt=dt)
        
        # Select protein
        protein = u.select_atoms("protein and name CA")  # C-alpha atoms
        
        if len(u.trajectory) < 3:
            raise ValueError(f"Insufficient frames for RMSF calculation: {len(u.trajectory)}")
        
        # Calculate RMSF
        R = RMSF(protein).run()
        rmsf_values = R.results.rmsf
        
        # Get residue IDs
        if original_resids is not None and len(original_resids) == len(protein.residues):
            resids = original_resids
        else:
            resids = protein.residues.resids
        
        # If binding pocket specified, filter to those residues
        if binding_pocket_resids is not None:
            # Convert binding_pocket_resids to integers for comparison
            pocket_resids_int = [int(r) for r in binding_pocket_resids]
            
            # Create mask for binding pocket residues
            mask = np.isin(resids, pocket_resids_int)
            
            filtered_resids = resids[mask]
            filtered_rmsf = rmsf_values[mask]
        else:
            filtered_resids = resids
            filtered_rmsf = rmsf_values
        
        # Define flexibility threshold (mean + 1 std)
        mean_rmsf = np.mean(filtered_rmsf)
        std_rmsf = np.std(filtered_rmsf)
        threshold = mean_rmsf + std_rmsf
        
        # Identify flexible residues
        flexible_mask = filtered_rmsf > threshold
        flexible_residues = filtered_resids[flexible_mask]
        flexible_rmsf = filtered_rmsf[flexible_mask]
        
        return {
            'resids': filtered_resids,
            'rmsf_values': filtered_rmsf,
            'flexible_residues': flexible_residues,
            'flexible_rmsf': flexible_rmsf,
            'threshold': threshold,
            'mean': mean_rmsf,
            'std': std_rmsf,
            'all_resids': resids,
            'all_rmsf': rmsf_values
        }
    
    finally:
        if u is not None:
            u.trajectory.close()
            del u


def plot_rmsf(rmsf_data, output_path, structure_name, pocket_only=False):
    """
    Create RMSF plot with flexible residues highlighted
    
    Parameters:
    -----------
    rmsf_data : dict
        Output from calculate_rmsf()
    output_path : str
        Path to save plot
    structure_name : str
        Name for plot title
    pocket_only : bool
        Whether this is pocket-only analysis
    """
    plt.figure(figsize=(12, 6))
    
    resids = rmsf_data['resids']
    rmsf_values = rmsf_data['rmsf_values']
    threshold = rmsf_data['threshold']
    
    # Plot RMSF
    plt.plot(resids, rmsf_values, 'o-', linewidth=1.5, markersize=4, label='RMSF', zorder=2)
    
    # Highlight flexible residues
    flexible_res = rmsf_data['flexible_residues']
    flexible_rmsf = rmsf_data['flexible_rmsf']
    plt.scatter(flexible_res, flexible_rmsf, color='red', s=50, 
               label=f'Flexible residues (n={len(flexible_res)})', zorder=5)
    
    # Add threshold line
    plt.axhline(y=threshold, color='r', linestyle='--', alpha=0.7,
               label=f'Flexibility threshold: {threshold:.3f} Å')
    
    # Add mean line
    plt.axhline(y=rmsf_data['mean'], color='g', linestyle='--', alpha=0.7,
               label=f'Mean: {rmsf_data["mean"]:.3f} Å')
    
    plt.xlabel('Residue Number', fontsize=12)
    plt.ylabel('RMSF (Å…)', fontsize=12)
    title_suffix = " (Binding Pocket)" if pocket_only else ""
    plt.title(f'Residue Flexibility Analysis: {structure_name}{title_suffix}', fontsize=14)
    plt.legend(fontsize=9, loc='best')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def calculate_relaxation_rmsd(trajectory_path, dt=5.0):
    """
    Calculate Cα RMSD of each frame versus the first frame, WITH rigid-body
    superposition, as a QC curve for how far the structure relaxes from the
    starting (predicted) model and whether it settles.

    Superposition (least-squares fit on Cα) is essential: the simulation runs in
    implicit solvent with no periodic box, so the protein is free to rotate and
    translate. Without alignment that global tumbling leaks into the RMSD and
    masquerades as conformational change. This curve is used for QC/plotting
    only — it no longer gates ensemble selection.

    Parameters
    ----------
    trajectory_path : str
        Path to the trajectory PDB file.
    dt : float
        Time between frames in ps.

    Returns
    -------
    dict: 'time' array and 'rmsd' array (Å, superposed vs frame 0).
    """
    u = None
    try:
        u = mda.Universe(trajectory_path, dt=dt)
        n_frames = len(u.trajectory)

        if n_frames < 2:
            raise ValueError(f"Trajectory has only {n_frames} frame(s), need at least 2 for relaxation RMSD")

        # RMSD vs first frame with superposition on Cα (select + groupselections empty)
        u.trajectory[0]
        rmsd_analysis = RMSD(u, u, select="protein and name CA", ref_frame=0).run()
        results = rmsd_analysis.results.rmsd  # columns: frame, time (ps), RMSD (Å)

        return {
            'time': results[:, 1],
            'rmsd': results[:, 2],
        }

    finally:
        if u is not None:
            u.trajectory.close()
            del u


def pocket_rmsd_matrix(trajectory_path, pocket_resids, original_resids, start_index):
    """
    Build a pairwise Cα RMSD matrix over binding-pocket residues for the frames
    from start_index onward, after rigid-body superposition of the whole protein.

    Aligning on the whole-protein Cα first, then measuring RMSD on only the pocket
    Cα atoms, captures pocket motion *relative to the body* rather than global
    drift. This matrix is the input to medoid-based ensemble selection.

    Residue identity is resolved the same way as calculate_rmsf: when
    original_resids is supplied it overrides the trajectory's internal (1..N)
    numbering so pocket_resids (in the reference numbering) select the right
    residues. If pocket_resids is falsy, all Cα atoms are used.

    Returns
    -------
    (frame_indices, matrix) : (list[int], np.ndarray)
        Global frame indices included and the symmetric NxN RMSD matrix (Å).
    """
    # Build the Universe without a `dt` kwarg: MDAnalysis forwards `dt` to the
    # in-memory MemoryReader (which already sets it), raising "got multiple values
    # for keyword argument 'dt'". Time is irrelevant here — only frame coordinates
    # matter — so AlignTraj(in_memory=True) does the in-memory transfer instead.
    u = mda.Universe(trajectory_path)
    try:
        # Superpose every frame onto frame 0 on protein Cα (removes tumbling/drift)
        align.AlignTraj(u, u, select="protein and name CA", ref_frame=0, in_memory=True).run()

        ca = u.select_atoms("protein and name CA")
        if original_resids is not None and len(original_resids) == len(ca.residues):
            resids = np.asarray(original_resids)
        else:
            resids = ca.residues.resids

        if pocket_resids:
            mask = np.isin(resids, [int(r) for r in pocket_resids])
            if not mask.any():
                mask = np.ones(len(ca), dtype=bool)  # fall back to all Cα
        else:
            mask = np.ones(len(ca), dtype=bool)

        n_frames = len(u.trajectory)
        frame_indices = list(range(start_index, n_frames))

        # Collect aligned pocket-Cα coordinates for each retained frame
        coords = []
        for idx in frame_indices:
            u.trajectory[idx]
            coords.append(ca.positions[mask].copy())
        coords = np.asarray(coords)

        n = len(coords)
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                diff = coords[i] - coords[j]
                rmsd = np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
                matrix[i, j] = matrix[j, i] = rmsd

        return frame_indices, matrix
    finally:
        u.trajectory.close()
        del u


def extract_ensemble_structures(trajectory_path, output_dir, pdb_name, frame_indices,
                                 original_resids=None,
                                 preserve_ligand=False, ligand_name=None, original_pdb=None,
                                 align_reference=None):
    """
    Write the given (medoid) trajectory frames as an ensemble of PDB structures.

    Before writing, every frame is rigid-body superposed (protein Cα) onto
    align_reference — the pre-MD input structure, which sits in the reference /
    docking-gridbox frame established by multiple_prot_align and shares the
    trajectory's topology (so the fit is exact). This puts the written structures in
    that frame rather than letting them carry the MD's rigid-body tumbling (in
    implicit solvent with no periodic box, OpenMM removes COM translation but not
    rotation). Otherwise each pocket drifts off the single reference-derived gridbox
    and the absolute-coordinate metrics (lig_RMSD / PLIF / PPS) are inflated. The
    internal conformational differences between medoids are preserved. Falls back to
    frame 0 if align_reference is missing or its atoms don't match.
    """
    u = None
    ensemble_files = []

    # Read ligand lines from original PDB if needed
    ligand_lines = []
    if preserve_ligand and ligand_name and original_pdb:
        with open(original_pdb, 'r') as f:
            original_lines = f.readlines()
        ligand_lines = [line for line in original_lines
                       if line.startswith('HETATM') and line[17:20].strip() == ligand_name]

    try:
        u = mda.Universe(trajectory_path)
        n_traj = len(u.trajectory)
        frame_indices = [idx for idx in frame_indices if 0 <= idx < n_traj]
        if not frame_indices:
            raise ValueError("No valid frame indices supplied for ensemble extraction")

        # Remove MD rigid-body drift: superpose every frame (protein Cα) onto the
        # reference frame. align_reference (the pre-MD input structure) is in the
        # multiple_prot_align / gridbox frame and shares the trajectory's topology, so
        # the fit is exact; frame 0 is the fallback when no reference is available.
        ref_universe = None
        if align_reference and os.path.exists(align_reference):
            try:
                ref_universe = mda.Universe(align_reference)
                align.AlignTraj(u, ref_universe, select="protein and name CA",
                                in_memory=True).run()
            except Exception as e:
                print(f"  Warning: could not align ensemble to reference structure "
                      f"({e}); falling back to frame 0")
                ref_universe = None
        if ref_universe is None:
            align.AlignTraj(u, u, select="protein and name CA", ref_frame=0,
                            in_memory=True).run()

        protein = u.select_atoms("protein")
        base_name = os.path.basename(pdb_name).replace(".pdb", "")
        
        for i, frame_idx in enumerate(frame_indices):
            u.trajectory[frame_idx]
            ensemble_filename = os.path.join(output_dir, f"{base_name}_ensemble_{i+1}.pdb")
            
            # Write structure with original residue numbering if provided
            with mda.Writer(ensemble_filename, protein.n_atoms) as W:
                if original_resids is not None:
                    if len(original_resids) == len(protein.residues):
                        protein.residues.resids = original_resids
                    else:
                        print(f"Warning: Original resid count mismatch for ensemble {i+1}. Using current numbering.")
                W.write(protein)
            
            # Append ligand if requested
            if preserve_ligand and ligand_lines:
                with open(ensemble_filename, 'r') as f:
                    protein_lines = f.readlines()
                
                # Insert ligand before TER/END
                insert_index = len(protein_lines)
                for idx, line in enumerate(protein_lines):
                    if line.startswith('TER') or line.startswith('END'):
                        insert_index = idx
                        break
                
                new_lines = protein_lines[:insert_index] + ligand_lines + protein_lines[insert_index:]
                
                with open(ensemble_filename, 'w') as f:
                    f.writelines(new_lines)
            
            ensemble_files.append(ensemble_filename)
        
        if preserve_ligand and ligand_lines:
            print(f"  Ligand {ligand_name} preserved in {len(ensemble_files)} ensemble structures")
        
        return ensemble_files
    
    finally:
        if u is not None:
            u.trajectory.close()
            del u


def run_mds(processed_filename, output_dir, mds_time=MDS_TIME_PS, original_resids=None,
            analyze_flexibility=False, binding_pocket_resids=None, structure_name=None,
            preserve_ligand=False, ligand_name=None, original_pdb=None):
    """
    Run Molecular Dynamics Simulation on a PDB file

    Runs a fixed-length relaxation simulation, then extracts an ensemble of
    N_ENSEMBLE cluster medoids from the post-burn-in frames (sampling within-basin
    conformational diversity) and profiles residue flexibility via RMSF.

    Parameters:
    -----------
    processed_filename : str
        Path to the processed PDB file
    output_dir : str
        Directory to save output files
    mds_time : int, optional
        Total simulation time in ps (defaults to the fixed MDS_TIME_PS)
    original_resids : np.array, optional
        Original residue IDs to preserve numbering
    analyze_flexibility : bool
        Whether to perform RMSF analysis for binding pocket flexibility
    binding_pocket_resids : list, optional
        List of binding pocket residue IDs for RMSF analysis and pocket-RMSD
        clustering
    structure_name : str, optional
        Name for output files (used in RMSF plots)
    preserve_ligand : bool
        Whether to preserve ligand from original structure
    ligand_name : str, optional
        Name of ligand to preserve
    original_pdb : str, optional
        Path to original PDB file containing ligand

    Returns:
    --------
    dict: Contains paths to output files and ensemble/relaxation statistics
        - 'mds_trajectory': Path to MDS trajectory file
        - 'burn_in_index' / 'burn_in_time': discarded relaxation transient
        - 'ensemble_structures': List of medoid ensemble structure paths
        - 'ensemble_frames': trajectory frame indices of the medoids
        - 'ensemble_diversity': mean pairwise pocket RMSD among medoids (Å)
        - 'rmsf_data': RMSF analysis results (if analyze_flexibility=True)
    """
    try:
        # Get the pdb_name
        pdb_name = os.path.basename(processed_filename).replace("_fixed.pdb", ".pdb")
        if structure_name is None:
            structure_name = os.path.basename(pdb_name).replace(".pdb", "")
        
        mds_output_name = os.path.join(output_dir, os.path.basename(pdb_name).replace(".pdb", "-mds.pdb"))
        
        # Perform Molecular Dynamics Simulation
        print("\nPerforming " + str(mds_time) + " ps MD simulation for equilibration...")
        
        # Load PDB file
        pdb = PDBFile(processed_filename)
        forcefield = ForceField('amber14-all.xml', 'implicit/gbn2.xml') 
        modeller = Modeller(pdb.topology, pdb.positions)
            
        system = forcefield.createSystem(modeller.topology,
                                    nonbondedMethod=CutoffNonPeriodic,
                                    nonbondedCutoff=1*nanometer,
                                    constraints=HBonds,
                                    hydrogenMass=1.5*amu)
        integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
        # Fixed seed → reproducible trajectories (subject to platform/thread determinism)
        integrator.setRandomNumberSeed(INTEGRATOR_SEED)
        simulation = Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)

        # Energy minimization
        print("Performing energy minimization...")
        simulation.minimizeEnergy(maxIterations=1000)

        # Calculate steps and reporting interval
        steps_per_ps = 250  # for 4 fs timestep
        report_interval_steps = 1250  # Save every 5 ps (improved from 10 ps)
        reporting_dt_ps = report_interval_steps / steps_per_ps  # Calculate actual time between frames
        total_steps = int(mds_time * steps_per_ps)

        # Setup reporters
        simulation.reporters.append(PDBReporter(mds_output_name, report_interval_steps))
        simulation.reporters.append(StateDataReporter(stdout, report_interval_steps, step=True,
                potentialEnergy=True, temperature=True))

        # Production run
        print(f"\nStarting production run ({mds_time} ps)...")
        print("Total steps:", total_steps)
        simulation.step(total_steps)

        # Clean up simulation objects
        del simulation
        del pdb

        # ------------------------------------------------------------------
        # Relaxation QC + ensemble selection (no equilibration gating)
        # ------------------------------------------------------------------
        print("\nAnalyzing relaxation (superposed Cα RMSD vs first frame)...\n")

        # Validate that trajectory was created and has frames
        if not os.path.exists(mds_output_name):
            raise RuntimeError(f"MD trajectory file not created: {mds_output_name}")

        u_check = mda.Universe(mds_output_name)
        n_frames = len(u_check.trajectory)
        u_check.trajectory.close()
        del u_check

        if n_frames < 3:
            raise RuntimeError(f"Insufficient frames in trajectory: {n_frames} (need at least 3)")

        print(f"Trajectory contains {n_frames} frames")

        # Relaxation curve — QC only (superposed RMSD vs frame 0)
        relax = calculate_relaxation_rmsd(mds_output_name, dt=reporting_dt_ps)
        time_array = relax['time']
        rmsd_values = relax['rmsd']

        # Burn-in: discard the leading relaxation transient before sampling
        burn_in_index = apply_burn_in(n_frames, BURN_IN_FRACTION)
        burn_in_time = time_array[burn_in_index] if burn_in_index < len(time_array) else time_array[-1]

        # Relaxation QC plot. The curve shows how far the structure moves from the
        # starting model; it is descriptive only — no attempt is made to infer
        # "equilibration", which is not achievable on this timescale.
        plt.figure(figsize=(10, 6))
        plt.plot(time_array, rmsd_values, marker='o', label='Cα RMSD vs frame 0 (superposed)', linewidth=1.5)
        plt.axvline(x=burn_in_time, color='r', linestyle='--',
                    label=f'Burn-in cutoff: {burn_in_time:.1f} ps ({int(BURN_IN_FRACTION*100)}%)')
        plt.xlabel("Time (ps)", fontsize=12)
        plt.ylabel(r'C$\alpha$ RMSD vs frame 0 ($\AA$)', fontsize=12)
        plt.title(f'Relaxation QC: {structure_name}', fontsize=14)
        plt.legend(fontsize=9)
        plt.grid(alpha=0.3)
        rmsd_plot_dir = os.path.join(output_dir, 'rmsd_plots')
        os.makedirs(rmsd_plot_dir, exist_ok=True)
        rmsd_plot_name = os.path.join(rmsd_plot_dir, f'{structure_name}_relaxation.png')
        plt.savefig(rmsd_plot_name, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"\nRelaxation QC for {structure_name}:")
        print(f"  Final Cα RMSD vs start: {rmsd_values[-1]:.3f} Å")
        print(f"  Burn-in discarded: first {burn_in_index} frame(s) (≈{burn_in_time:.1f} ps)")

        # Pocket-RMSD matrix over post-burn-in frames → medoid ensemble selection
        frame_pool, dmatrix = pocket_rmsd_matrix(
            mds_output_name, binding_pocket_resids, original_resids, burn_in_index
        )
        medoid_local = select_medoids(dmatrix, N_ENSEMBLE)
        medoid_frames = [frame_pool[i] for i in medoid_local]
        diversity = intra_ensemble_diversity(dmatrix, medoid_local)
        
        # Write the medoid ensemble. These cluster medoids replace both the old
        # single representative and the linspace-over-plateau ensemble: each test
        # structure is now represented by N_ENSEMBLE diverse, relaxed conformations.
        print(f"\nExtracting {len(medoid_frames)}-structure ensemble (cluster medoids)...")
        ensemble_files = extract_ensemble_structures(
            mds_output_name,
            output_dir,
            pdb_name,
            medoid_frames,
            original_resids=original_resids,
            preserve_ligand=preserve_ligand,
            ligand_name=ligand_name,
            original_pdb=original_pdb,
            align_reference=processed_filename,
        )
        print(f"  Extracted {len(ensemble_files)} ensemble structure(s) from frames {medoid_frames}")
        print(f"  Intra-ensemble diversity (mean pairwise pocket RMSD): {diversity:.3f} Å")
        if diversity < LOW_DIVERSITY_RMSD:
            print(f"  ⚠ Low diversity (< {LOW_DIVERSITY_RMSD} Å): pocket is effectively rigid; "
                  f"ensemble members are near-duplicates")

        # Perform RMSF analysis if requested
        rmsf_data = None
        if analyze_flexibility:
            print("\nPerforming RMSF analysis for binding pocket flexibility...")
            
            try:
                rmsf_data = calculate_rmsf(
                    mds_output_name,
                    binding_pocket_resids=binding_pocket_resids,
                    dt=reporting_dt_ps,
                    original_resids=original_resids
                )
                
                # Create RMSF plots directory
                rmsf_plot_dir = os.path.join(output_dir, 'rmsf_plots')
                os.makedirs(rmsf_plot_dir, exist_ok=True)
                
                # Plot binding pocket RMSF
                if binding_pocket_resids:
                    pocket_plot_path = os.path.join(rmsf_plot_dir, 
                                                   f'{structure_name}_rmsf_pocket.png')
                    plot_rmsf(rmsf_data, pocket_plot_path, structure_name, pocket_only=True)
                    
                    print(f"  Binding pocket residues analyzed: {len(rmsf_data['resids'])}")
                    print(f"  Flexible residues identified: {len(rmsf_data['flexible_residues'])}")
                    print(f"  Flexible residues: {sorted(rmsf_data['flexible_residues'])}")
                
                # Also create full protein plot
                full_plot_path = os.path.join(rmsf_plot_dir, 
                                             f'{structure_name}_rmsf_full.png')
                full_rmsf_data = calculate_rmsf(mds_output_name, 
                                               binding_pocket_resids=None,
                                               dt=reporting_dt_ps,
                                               original_resids=original_resids)
                plot_rmsf(full_rmsf_data, full_plot_path, structure_name, pocket_only=False)
                
            except Exception as e:
                print(f"Warning: RMSF analysis failed: {e}")
                import traceback
                traceback.print_exc()

        # Clean up processed file after successful completion
        if os.path.exists(processed_filename):
            try:
                os.remove(processed_filename)
            except OSError as e:
                print(f"Warning: Could not remove {processed_filename}: {e}")

        # Return comprehensive results
        return {
            'structure_name': structure_name,
            'mds_trajectory': mds_output_name,
            'burn_in_index': burn_in_index,
            'burn_in_time': burn_in_time,
            'ensemble_structures': ensemble_files,
            'ensemble_frames': medoid_frames,
            'ensemble_diversity': diversity,
            'rmsf_data': rmsf_data
        }

    except Exception as e:
        # Clean up on error
        print(f"Error during MDS: {str(e)}")
        
        if os.path.exists(processed_filename):
            try:
                os.remove(processed_filename)
            except OSError as oe:
                print(f"Warning: Failed to clean up {processed_filename}: {oe}")
        
        raise


def prep_receptors():
    """
    Process all PDB files in a specified directory
    
    Interactive workflow for preparing protein structures for docking:
    1. Fix structures with PDBFixer
    2. Optionally run MD simulations for equilibration
    3. Extract equilibrated structures with frame-to-frame RMSD analysis
    4. For reference structure: extract ensemble for conformational sampling
    5. For test structures: analyze binding pocket flexibility via RMSF
    """
    # Set working directory
    pdb_dir = input("Enter the path to the directory containing the aligned PDB files: ")
    while not os.path.exists(pdb_dir):
        pdb_dir = input(
            "That path does not appear to exist.\nPlease enter the path to "
            "the directory containing the PDB files: "
        )
    pdb_dir = os.path.abspath(pdb_dir)
    os.chdir(pdb_dir)

    # Create output directory
    output_dir = 'prepared_structures'
    os.makedirs(output_dir, exist_ok=True)
    
    # Handle reference structure - check for files starting with "ref_"
    ref_name = None
    ref_ligand = None
    binding_pocket_resids = None
    
    # Search for files starting with "ref_"
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    ref_candidates = [f for f in pdb_files if f.startswith('ref_')]
    
    if ref_candidates:
        print(f"\nFound {len(ref_candidates)} file(s) starting with 'ref_':")
        for candidate in ref_candidates:
            if prompt_yes_no(f"Is {candidate} the reference PDB file? (y/n): "):
                ref_name = candidate
                break

        if not ref_name:
            print("No reference structure selected from candidates.")
            if prompt_yes_no("Would you like to manually specify a reference structure? (y/n): "):
                ref_name = input("Enter the name of the reference structure: ")
                if not ref_name.endswith(".pdb"):
                    ref_name += ".pdb"
                if ref_name not in pdb_files:
                    print(f"Reference structure {ref_name} not found in the directory.")
                    ref_name = None
    else:
        # No ref_ files found, ask user
        if prompt_yes_no("Is there a reference structure (with ligand bound) in the specified PDB directory? (y/n): "):
            ref_name = input("Enter the name of the reference structure: ")
            if not ref_name.endswith(".pdb"):
                ref_name += ".pdb"
            if ref_name not in pdb_files:
                print(f"Reference structure {ref_name} not found in the directory.")
                ref_name = None
    
    # If we have a reference, get ligand info and identify binding pocket
    if ref_name:
        ref_path = os.path.join(pdb_dir, ref_name)

        # Detect candidate ligand residue names from the reference HETATM records
        # and let the user pick one, rather than typing it blind.
        candidates = list_hetatm_ligands(ref_path)
        if candidates:
            print(f"\nHETATM residues found in {ref_name} (largest first):")
            for i, (resname, n_atoms) in enumerate(candidates, 1):
                print(f"  [{i}] {resname}  ({n_atoms} atoms)")
            print("  [m] Enter a residue name manually")
            while True:
                choice = input("Select the reference ligand (number, or 'm'): ").strip().lower()
                if choice == "m":
                    ref_ligand = input("Enter the 3-letter ligand code: ").strip()
                    break
                if choice.isdigit() and 1 <= int(choice) <= len(candidates):
                    ref_ligand = candidates[int(choice) - 1][0]
                    break
                print("  Invalid selection, please try again.")
        else:
            print(f"\nNo HETATM records found in {ref_name}.")
            ref_ligand = input("Enter the 3-letter ligand code in the reference structure: ").strip()
        print(f"  Using reference ligand: {ref_ligand}")

        # Identify binding pocket residues
        print(f"\nIdentifying binding pocket residues in {ref_name}...")
        binding_pocket_resids = get_binding_pocket_residues(ref_path, ref_ligand)
        
        if binding_pocket_resids:
            print(f"  Found {len(binding_pocket_resids)} binding pocket residues: {binding_pocket_resids[:10]}{'...' if len(binding_pocket_resids) > 10 else ''}")
        else:
            print(f"  Warning: No binding pocket residues identified")

    # Get all PDB files in directory
    all_pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]

    # Separate reference from test structures
    test_pdb_files = [f for f in all_pdb_files if f != ref_name] if ref_name else all_pdb_files
    
    print(f"\nFound {len(test_pdb_files)} test structure(s) and {1 if ref_name else 0} reference structure")
    
    perform_mds = prompt_yes_no("Perform a short (200 ps) Molecular Dynamics Simulation for structure relaxation? (y/n): ")

    if not perform_mds:
        # Just run PDBFixer
        print("\nRunning PDBFixer on all PDB files...")
        results = []
        
        # Process reference if present
        if ref_name:
            full_path = os.path.join(pdb_dir, ref_name)
            try:
                fixed_file, _ = run_pdbfixer(full_path, output_dir, 
                                            keep_ligand=True, 
                                            ligand_name=ref_ligand)
                
                # Determine final reference name - add ref_ prefix only if not already present
                if ref_name.startswith('ref_'):
                    ref_output = os.path.join(output_dir, ref_name)
                else:
                    ref_output = os.path.join(output_dir, f"ref_{ref_name}")
                
                # Copy to output with appropriate name
                shutil.copy(fixed_file, ref_output)
                os.remove(fixed_file)
                results.append({'fixed_structure': ref_output, 'is_reference': True})
            except Exception as e:
                print(f"Failed to process reference {ref_name}. Error: {e}")
        
        # Process test structures
        for pdb_file in test_pdb_files:
            full_path = os.path.join(pdb_dir, pdb_file)
            try:
                fixed_file, _ = run_pdbfixer(full_path, output_dir)
                results.append({'fixed_structure': fixed_file, 'is_reference': False})
            except Exception as e:
                print(f"Failed to process {pdb_file} with PDBFixer. Error: {e}")
        
        print("\nAll PDB files have been processed with PDBFixer.")
        return results
    
    # MD simulation parameters.
    # The simulation length, ensemble size, and flexibility analysis are no longer
    # user-tunable: the protocol is standardized (fixed MDS_TIME_PS, N_ENSEMBLE medoid
    # structures, flexibility profiling always on) for reproducibility. The only
    # remaining choice is whether to keep the (large) trajectory files.
    save_trajectories = prompt_yes_no("Do you want to save the MDS trajectory files? (y/n): ")

    # Flexibility profiling is always performed when a binding pocket is available.
    analyze_test_flexibility = bool(binding_pocket_resids)

    print(f"\nStandardized MD protocol: {MDS_TIME_PS} ps, {N_ENSEMBLE} medoid ensemble "
          f"structures per receptor, seed={INTEGRATOR_SEED}.")

    results = []
    
    # Process REFERENCE structure first (if present)
    if ref_name:
        ref_path = os.path.join(pdb_dir, ref_name)
        try:
            print(f"\n{'='*60}")
            print(f"Processing REFERENCE structure: {ref_name}")
            print(f"  (Using crystal structure - no MD simulation)")
            print('='*60)
            
            # Process reference structure with ligand preservation (no MD)
            processed_filename, original_resids = run_pdbfixer(ref_path, output_dir, 
                                                               keep_ligand=True, 
                                                               ligand_name=ref_ligand)
            
            # Determine final reference name - add ref_ prefix only if not already present
            if ref_name.startswith('ref_'):
                final_ref_name = os.path.join(output_dir, ref_name)
            else:
                final_ref_name = os.path.join(output_dir, f"ref_{ref_name}")
            
            # Move the fixed file to final location
            shutil.move(processed_filename, final_ref_name)
            
            print(f"\n  Reference structure prepared: {os.path.basename(final_ref_name)}")
            print(f"  Contains ligand: {ref_ligand}")
            print(f"  Source: Crystal structure (no MD equilibration)")
            
            # Create result dictionary
            result = {
                'is_reference': True,
                'ligand': ref_ligand,
                'binding_pocket_resids': binding_pocket_resids,
                'reference_structure': final_ref_name,
                'method': 'crystal_structure'
            }
            results.append(result)
        
        except Exception as e:
            print(f"Failed to process reference {ref_name}. Error: {e}")
            import traceback
            traceback.print_exc()
    
    # Process TEST structures
    print(f"\n{'='*60}")
    print(f"Processing {len(test_pdb_files)} TEST structures")
    print('='*60)
    
    for pdb_file in test_pdb_files:
        full_path = os.path.join(pdb_dir, pdb_file)
        try:
            print(f"\n{'='*60}")
            print(f"Processing: {pdb_file}")
            print('='*60)
            
            processed_filename, original_resids = run_pdbfixer(full_path, output_dir)
            result = run_mds(
                processed_filename,
                output_dir,
                MDS_TIME_PS,
                original_resids=original_resids,
                analyze_flexibility=analyze_test_flexibility,
                binding_pocket_resids=binding_pocket_resids,
                structure_name=pdb_file.replace('.pdb', '')
            )
            
            result['is_reference'] = False
            results.append(result)
            
            # Clean up trajectory if not saving
            if not save_trajectories and os.path.exists(result['mds_trajectory']):
                try:
                    os.remove(result['mds_trajectory'])
                    print(f"Removed trajectory file: {os.path.basename(result['mds_trajectory'])}")
                except OSError as e:
                    print(f"Warning: Could not remove trajectory {result['mds_trajectory']}: {e}")
        
        except Exception as e:
            print(f"Failed to process {pdb_file}. Error: {e}")
            import traceback
            traceback.print_exc()

    # Print summary
    print("\n" + "="*60)
    print("Processing Complete")
    print("="*60)
    
    # Use builtin sum explicitly to avoid OpenMM override conflict
    ref_count = __builtins__.sum(1 for r in results if r.get('is_reference', False))
    test_count = len(results) - ref_count
    
    print(f"\nSuccessfully processed:")
    print(f"  Reference structures: {ref_count}")
    print(f"  Test structures: {test_count}")
    print(f"  Total: {len(results)}")
    print(f"\nOutput directory: {os.path.abspath(output_dir)}")
    
    # Summary of reference structure
    if ref_count > 0:
        ref_result = [r for r in results if r.get('is_reference', False)][0]
        print(f"\nReference structure: {os.path.basename(ref_result['reference_structure'])}")
        print(f"  Method: Crystal structure (no MD)")
        print(f"  Ligand: {ref_result.get('ligand', 'N/A')}")
        print(f"  Binding pocket residues: {len(ref_result.get('binding_pocket_resids', []))}")
    
    # Summary of flexibility analysis
    if analyze_test_flexibility:
        flex_analyzed = [r for r in results if r.get('rmsf_data') is not None]
        print(f"\nFlexibility analysis completed for {len(flex_analyzed)} test structures")
        
        for result in flex_analyzed:
            if result.get('rmsf_data'):
                rmsf_data = result['rmsf_data']
                structure = result['structure_name']
                print(f"  {structure}: {len(rmsf_data['flexible_residues'])} flexible residues")
        
        # Create consolidated flexible residues summary file
        details_dir = os.path.join(pdb_dir, "results")
        save_flexible_residues_summary(details_dir, results)

    # Persist project-level facts so downstream scripts can skip re-prompting
    config = load_config(pdb_dir)
    config["project_dir"] = pdb_dir
    if ref_name:
        config["reference_pdb"] = ref_name
    if ref_ligand:
        config["ligand_resname"] = ref_ligand
    save_config(pdb_dir, config)
    print(f"\nProject config saved to: {os.path.join(pdb_dir, CONFIG_FILENAME)}")

    return results


if __name__ == "__main__":
    prep_receptors()