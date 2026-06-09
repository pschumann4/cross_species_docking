"""
mds_structure_prep.py
=============================
Prepares the test and reference structure for molecular dynamics simulation and flexibility analysis. 
It performs the following steps:
1. Runs PDBFixer to clean up the structure, add missing atoms, and ensure it's suitable for MD.
2. Identifies binding pocket residues based on proximity to the ligand (in reference structure).
3. Runs a short MD simulation to equilibrate the structure and generate a trajectory.
4. Analyzes the trajectory to calculate RMSF (Root Mean Square Fluctuation) for each residue, 
   identifying which residues are flexible based on a defined threshold.
5. Optionally extracts multiple structures from the equilibrated plateau region to create an ensemble.

Outputs
-------
- Fixed PDB files for test and reference structures (protein only for MD, with ligand for reference)
- RMSF plots showing residue flexibility with flexible residues highlighted
- A summary file listing flexible residues for each structure
- Extracted ensemble structures from the plateau region (if enabled)
"""

import os
import sys
import shutil
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import euclidean3d, centroid, load_config, save_config, CONFIG_FILENAME
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.rms import RMSF
import ruptures as rpt
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from pdbfixer import PDBFixer
from openmm.app import PDBFile


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
            structure_name = os.path.basename(result['equilibration_structure']).replace('.pdb', '')
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


def calculate_frame_to_frame_rmsd(trajectory_path, dt=5.0):
    """
    Calculate frame-to-frame RMSD to detect when structure stops changing
    
    Parameters:
    -----------
    trajectory_path : str
        Path to the trajectory PDB file
    dt : float
        Time between frames in ps (default: 5.0)
    
    Returns:
    --------
    dict: Contains 'time' array and 'rmsd' array (frame-to-frame RMSD values)
    """
    u = None
    try:
        u = mda.Universe(trajectory_path, dt=dt)
        n_frames = len(u.trajectory)
        
        if n_frames < 2:
            raise ValueError(f"Trajectory has only {n_frames} frame(s), need at least 2 for frame-to-frame RMSD")
        
        protein = u.select_atoms("protein and name CA")  # Use C-alpha atoms for efficiency
        
        # Initialize arrays
        frame_to_frame_rmsd = []
        time_points = []
        
        # Set reference to first frame
        u.trajectory[0]
        ref_positions = protein.positions.copy()
        
        # Calculate RMSD between consecutive frames
        for ts in u.trajectory[1:]:
            current_positions = protein.positions.copy()
            
            # Calculate RMSD between current and previous frame
            diff = current_positions - ref_positions
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            
            frame_to_frame_rmsd.append(rmsd)
            time_points.append(ts.time)
            
            # Update reference to current frame
            ref_positions = current_positions
        
        return {
            'time': np.array(time_points),
            'rmsd': np.array(frame_to_frame_rmsd)
        }
    
    finally:
        if u is not None:
            u.trajectory.close()
            del u


def estimate_plateau_point(rmsd_values, time, min_plateau_duration_ps=20.0, pen_factor=2.0):
    """
    Estimate the equilibration point using changepoint detection
    
    Uses PELT algorithm to detect when RMSD stabilizes. Includes validation
    to ensure the plateau is sustained for a minimum duration.
    
    Parameters:
    -----------
    rmsd_values : np.array
        Array of RMSD values
    time : np.array
        Array of time points corresponding to RMSD values
    min_plateau_duration_ps : float
        Minimum duration (in ps) that a plateau must be sustained to be valid
    pen_factor : float
        Penalty factor for PELT algorithm (higher = fewer changepoints)
    
    Returns:
    --------
    dict: Plateau statistics including start time, mean, std, and representative frame
    """
    if len(rmsd_values) < 3:
        raise ValueError(f"Need at least 3 data points for plateau detection, got {len(rmsd_values)}")
    
    # Adaptive penalty based on data characteristics
    data_range = np.ptp(rmsd_values)  # Peak-to-peak range
    penalty = pen_factor * data_range
    
    # Plateau detection using PELT algorithm
    try:
        algo = rpt.Pelt(model="rbf").fit(rmsd_values)
        changepoints = algo.predict(pen=penalty)
    except Exception as e:
        print(f"Warning: PELT algorithm failed with pen={penalty:.3f}, using fallback method")
        # Fallback: use the point where RMSD drops below mean + 0.5*std
        mean_rmsd = np.mean(rmsd_values)
        std_rmsd = np.std(rmsd_values)
        threshold = mean_rmsd + 0.5 * std_rmsd
        
        below_threshold = rmsd_values < threshold
        if np.any(below_threshold):
            plateau_start_index = np.argmax(below_threshold)
        else:
            # Last resort: use midpoint
            plateau_start_index = len(rmsd_values) // 2
        changepoints = [plateau_start_index]
    
    # Use first changepoint as potential plateau start
    plateau_start_index = changepoints[0] if changepoints else len(rmsd_values) // 2
    
    # Ensure we don't start beyond the data
    if plateau_start_index >= len(rmsd_values):
        plateau_start_index = len(rmsd_values) // 2
        print(f"Warning: Changepoint beyond data range, using midpoint at index {plateau_start_index}")
    
    # Validate minimum plateau duration
    plateau_duration = time[-1] - time[plateau_start_index]
    
    if plateau_duration < min_plateau_duration_ps:
        print(f"Warning: Detected plateau duration ({plateau_duration:.1f} ps) is less than "
              f"minimum required ({min_plateau_duration_ps:.1f} ps)")
        # Try to find an earlier stable region
        target_index = len(rmsd_values) - int(min_plateau_duration_ps / (time[1] - time[0]))
        if target_index < 0:
            print("Warning: Insufficient simulation time for minimum plateau duration")
            target_index = 0
        plateau_start_index = max(0, target_index)
    
    # Calculate plateau statistics
    plateau_values = rmsd_values[plateau_start_index:]
    if len(plateau_values) == 0:
        raise ValueError("No data points in plateau region")
    
    plateau_average = np.mean(plateau_values)
    plateau_std = np.std(plateau_values)
    
    # Find the frame closest to the plateau mean
    differences = np.abs(plateau_values - plateau_average)
    closest_to_mean_local_idx = np.argmin(differences)
    closest_to_mean_index = plateau_start_index + closest_to_mean_local_idx

    return {
        'start_index': plateau_start_index,
        'start_time': time[plateau_start_index],
        'plateau_average': plateau_average,
        'plateau_std': plateau_std,
        'mean_representative_index': closest_to_mean_index,
        'mean_representative_time': time[closest_to_mean_index]
    }


def extract_ensemble_structures(trajectory_path, output_dir, pdb_name, plateau_start_index, 
                                 plateau_end_index=None, n_structures=5, original_resids=None,
                                 preserve_ligand=False, ligand_name=None, original_pdb=None):
    """
    Extract multiple structures from the equilibrated plateau region
    
    This provides an ensemble of structures capturing conformational diversity,
    which can compensate for any inaccuracies in the structural prediction and 
    provide better sampling for docking.
    
    Parameters:
    -----------
    trajectory_path : str
        Path to trajectory file
    output_dir : str
        Output directory
    pdb_name : str
        Base name for output files
    plateau_start_index : int
        First frame of plateau
    plateau_end_index : int, optional
        Last frame of plateau (default: end of trajectory)
    n_structures : int
        Number of structures to extract from plateau
    original_resids : np.array, optional
        Original residue IDs to preserve numbering
    preserve_ligand : bool
        Whether to append ligand from original structure
    ligand_name : str, optional
        Name of ligand to preserve
    original_pdb : str, optional
        Path to original PDB file containing ligand
    
    Returns:
    --------
    list: Paths to extracted PDB files
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
        
        if plateau_end_index is None:
            plateau_end_index = len(u.trajectory) - 1
        
        # Validate indices
        plateau_end_index = min(plateau_end_index, len(u.trajectory) - 1)
        
        if plateau_start_index >= len(u.trajectory):
            raise ValueError(f"Plateau start index {plateau_start_index} exceeds trajectory length {len(u.trajectory)}")
        
        # Select frames uniformly from plateau
        plateau_length = plateau_end_index - plateau_start_index + 1
        
        if plateau_length < n_structures:
            print(f"Warning: Plateau has only {plateau_length} frames, extracting all")
            frame_indices = list(range(plateau_start_index, plateau_end_index + 1))
        else:
            frame_indices = np.linspace(plateau_start_index, plateau_end_index, 
                                       n_structures, dtype=int)
        
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


def run_mds(processed_filename, output_dir, mds_time=None, extract_ensemble=False, 
            n_ensemble=1, min_plateau_duration=20.0, original_resids=None,
            analyze_flexibility=False, binding_pocket_resids=None, structure_name=None,
            preserve_ligand=False, ligand_name=None, original_pdb=None):
    """
    Run Molecular Dynamics Simulation on a PDB file

    Parameters:
    -----------
    processed_filename : str
        Path to the processed PDB file
    output_dir : str
        Directory to save output files
    mds_time : int, optional
        Total simulation time in ps
    extract_ensemble : bool
        Whether to extract multiple structures from plateau (default: False)
    n_ensemble : int
        Number of structures to extract if extract_ensemble=True (default: 1)
    min_plateau_duration : float
        Minimum plateau duration in ps for validation (default: 20.0)
    original_resids : np.array, optional
        Original residue IDs to preserve numbering
    analyze_flexibility : bool
        Whether to perform RMSF analysis for binding pocket flexibility
    binding_pocket_resids : list, optional
        List of binding pocket residue IDs for RMSF analysis
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
    dict: Contains paths to output files and equilibration statistics
        - 'mds_trajectory': Path to MDS trajectory file
        - 'equilibration_time': Time of equilibration in ps
        - 'equilibration_structure': Path to single equilibrated structure
        - 'ensemble_structures': List of ensemble structure paths (if extract_ensemble=True)
        - 'plateau_stats': Dictionary of plateau statistics
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

        # RMSD Analysis using frame-to-frame approach
        print("\nCalculating frame-to-frame RMSD to detect equilibration...\n")
        
        # Validate that trajectory was created and has frames
        if not os.path.exists(mds_output_name):
            raise RuntimeError(f"MD trajectory file not created: {mds_output_name}")
        
        # Quick check of frame count
        u_check = mda.Universe(mds_output_name)
        n_frames = len(u_check.trajectory)
        u_check.trajectory.close()
        del u_check
        
        if n_frames < 3:
            raise RuntimeError(f"Insufficient frames in trajectory: {n_frames} (need at least 3)")
        
        print(f"Trajectory contains {n_frames} frames")
        
        # Calculate frame-to-frame RMSD
        rmsd_results = calculate_frame_to_frame_rmsd(mds_output_name, dt=reporting_dt_ps)
        
        time_array = rmsd_results['time']
        rmsd_values = rmsd_results['rmsd']

        # Estimate plateau point with validation
        plateau_results = estimate_plateau_point(rmsd_values, time_array, 
                                                 min_plateau_duration_ps=min_plateau_duration,
                                                 pen_factor=2.0)
        
        plateau_index = plateau_results['mean_representative_index']
        plateau_time = plateau_results['start_time']

        # Plot RMSD with enhanced visualization
        plt.figure(figsize=(10, 6))
        plt.plot(time_array, rmsd_values, marker='o', label='Frame-to-Frame RMSD', linewidth=1.5)

        # Add plateau indicators
        plt.axvline(x=plateau_time, color='r', linestyle='--', 
                   label=f'Equilibration Start: {plateau_time:.1f} ps')
        plt.axvline(x=plateau_results['mean_representative_time'], color='b', 
                   linestyle='--', label=f'Representative Frame: {plateau_results["mean_representative_time"]:.1f} ps')
        plt.axhline(y=plateau_results['plateau_average'], 
                color='g', linestyle='--', 
                label=f'Plateau Mean: {plateau_results["plateau_average"]:.3f} Å')

        # Add plateau region shading
        plt.fill_between(time_array[plateau_results['start_index']:],
                        plateau_results['plateau_average'] - plateau_results['plateau_std'],
                        plateau_results['plateau_average'] + plateau_results['plateau_std'],
                        color='g', alpha=0.2,
                        label=f'Std Dev: ±{plateau_results["plateau_std"]:.3f} Å')

        plt.xlabel("Time (ps)", fontsize=12)
        plt.ylabel(r'Frame-to-Frame RMSD ($\AA$)', fontsize=12)
        plt.title(f'Equilibration Analysis: {structure_name}', fontsize=14)
        plt.legend(fontsize=9)
        plt.grid(alpha=0.3)
        
        # Save plot
        rmsd_plot_dir = os.path.join(output_dir, 'rmsd_plots')
        os.makedirs(rmsd_plot_dir, exist_ok=True)
        rmsd_plot_name = os.path.join(rmsd_plot_dir, f'{structure_name}_equilibration.png')
        plt.savefig(rmsd_plot_name, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"\nEquilibration analysis for {structure_name}:")
        print(f"  Representative frame time: {plateau_results['mean_representative_time']:.2f} ps")
        print(f"  Plateau average RMSD: {plateau_results['plateau_average']:.3f} Å")
        print(f"  Plateau std dev: {plateau_results['plateau_std']:.3f} Å")
        print(f"  Plateau duration: {time_array[-1] - plateau_time:.2f} ps")
        
        # Validate plateau_index is within bounds
        if plateau_index >= n_frames:
            print(f"Warning: Calculated plateau index {plateau_index} exceeds frame count {n_frames}")
            plateau_index = n_frames - 1
            print(f"Using last frame (index {plateau_index}) instead")
        
        # Extract single representative structure
        equilibration_pdb_name = os.path.join(output_dir, os.path.basename(pdb_name))
        u_eq = None
        try:
            u_eq = mda.Universe(mds_output_name)
            u_eq.trajectory[plateau_index]
            protein = u_eq.select_atoms("protein")
            
            # Use PDB writer with preserve option to keep original residue numbering
            with mda.Writer(equilibration_pdb_name, protein.n_atoms) as W:
                # If original resids provided, restore them before writing
                if original_resids is not None:
                    if len(original_resids) == len(protein.residues):
                        protein.residues.resids = original_resids
                    else:
                        print(f"Warning: Original resid count ({len(original_resids)}) doesn't match "
                              f"current residue count ({len(protein.residues)}). Using current numbering.")
                W.write(protein)
            
            # Append ligand if requested
            if preserve_ligand and ligand_name and original_pdb:
                with open(original_pdb, 'r') as f:
                    original_lines = f.readlines()
                ligand_lines = [line for line in original_lines 
                               if line.startswith('HETATM') and line[17:20].strip() == ligand_name]
                
                if ligand_lines:
                    with open(equilibration_pdb_name, 'r') as f:
                        protein_lines = f.readlines()
                    
                    # Insert ligand before TER/END
                    insert_index = len(protein_lines)
                    for idx, line in enumerate(protein_lines):
                        if line.startswith('TER') or line.startswith('END'):
                            insert_index = idx
                            break
                    
                    new_lines = protein_lines[:insert_index] + ligand_lines + protein_lines[insert_index:]
                    
                    with open(equilibration_pdb_name, 'w') as f:
                        f.writelines(new_lines)
                    
                    print(f"  Ligand {ligand_name} preserved in equilibrated structure")
            
            print(f"  Extracted representative structure: {os.path.basename(equilibration_pdb_name)}")
        finally:
            if u_eq is not None:
                u_eq.trajectory.close()
                del u_eq
        
        # Extract ensemble if requested
        ensemble_files = []
        if extract_ensemble and n_ensemble > 1:
            print(f"\nExtracting ensemble of {n_ensemble} structures from plateau...")
            ensemble_files = extract_ensemble_structures(
                mds_output_name, 
                output_dir, 
                pdb_name,
                plateau_results['start_index'],
                n_structures=n_ensemble,
                original_resids=original_resids,
                preserve_ligand=preserve_ligand,
                ligand_name=ligand_name,
                original_pdb=original_pdb
            )
            print(f"  Extracted {len(ensemble_files)} ensemble structures")

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
            'mds_trajectory': mds_output_name,
            'equilibration_time': plateau_time,
            'equilibration_structure': equilibration_pdb_name,
            'ensemble_structures': ensemble_files,
            'plateau_stats': plateau_results,
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
    pdb_dir = input("Enter the path to the directory containing the PDB files: ")
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
            confirm = input(f"Is {candidate} the reference PDB file? (y/n): ").lower().strip()
            if confirm == 'y':
                ref_name = candidate
                break
        
        if not ref_name:
            print("No reference structure selected from candidates.")
            manual_ref = input("Would you like to manually specify a reference structure? (y/n): ").lower().strip()
            if manual_ref == 'y':
                ref_name = input("Enter the name of the reference structure: ")
                if not ref_name.endswith(".pdb"):
                    ref_name += ".pdb"
                if ref_name not in pdb_files:
                    print(f"Reference structure {ref_name} not found in the directory.")
                    ref_name = None
    else:
        # No ref_ files found, ask user
        has_ref = input("Is there a reference structure (with ligand bound) in the specified PDB directory? (y/n): ").lower().strip()
        if has_ref == 'y':
            ref_name = input("Enter the name of the reference structure: ")
            if not ref_name.endswith(".pdb"):
                ref_name += ".pdb"
            if ref_name not in pdb_files:
                print(f"Reference structure {ref_name} not found in the directory.")
                ref_name = None
    
    # If we have a reference, get ligand info and identify binding pocket
    if ref_name:
        ref_path = os.path.join(pdb_dir, ref_name)
        
        # Ask for ligand name for binding pocket identification
        ref_ligand = input("Enter the 3-letter ligand code in the reference structure: ").strip()
        
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
    
    perform_mds = input("Do you want to perform Molecular Dynamics Simulation? (y/n): ").lower().strip() == 'y'

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
    
    # MD simulation parameters
    save_trajectories = input(f"Do you want to save the MDS trajectory files? (y/n): ").lower().strip() == 'y'
    
    mds_time = input("Enter the total simulation time in ps (suggest 100-500 ps): ")
    while not mds_time.isdigit() or int(mds_time) < 100 or int(mds_time) > 1000:
        mds_time = input("Please enter a valid time (100-1000 ps): ")
    mds_time = int(mds_time)
    
    # Test structure flexibility analysis
    analyze_test_flexibility = False
    if binding_pocket_resids:
        analyze_test_flexibility = input("Perform RMSF flexibility analysis on test structures' binding pockets? (y/n): ").lower().strip() == 'y'

    # Test structure ensemble extraction
    extract_test_ensemble = input("Do you want to save ensemble structures from test trajectories? (y/n): ").lower().strip() == 'y'
    n_test_ensemble = 1
    if extract_test_ensemble:
        n_ensemble_input = input("How many ensemble structures to extract per test structure? (suggest 3-5): ")
        while not n_ensemble_input.isdigit() or int(n_ensemble_input) < 2 or int(n_ensemble_input) > 10:
            n_ensemble_input = input("Please enter a valid number (2-10): ")
        n_test_ensemble = int(n_ensemble_input)

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
                mds_time,
                extract_ensemble=extract_test_ensemble,
                n_ensemble=n_test_ensemble,
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
                structure = os.path.basename(result['equilibration_structure']).replace('.pdb', '')
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