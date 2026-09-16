"""
Integrated AutoDock Vina workflow: batch docking → rank-1 pose extraction → PDB model
generation (obabel PDBQT→PDB conversion with reference-based HETATM ordering).
"""

import os
import sys
import csv
import subprocess
import glob
import shutil
import time
from datetime import datetime
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import check_tools, parse_hetatm_coords, rmsd_hungarian, prompt_yes_no

# vina_split numbers split poses consecutively from 1; this is the upper bound
# on how many modes to scan for / clean up (Vina's num_modes never approaches it).
MAX_VINA_MODES = 100

# Residue sequence number stamped onto the docked ligand in generated PDB models
# so the ligand is unambiguously distinct from protein residues. Must stay 4
# characters wide to fit PDB columns 23-26.
LIGAND_RESNUM = "9999"


# ============================================================================
# VINA BATCH PROCESSING FUNCTIONS
# ============================================================================

def run_vina(verbose=True):
    """
    Run AutoDock Vina for each configuration file in the current directory.
    
    Parameters:
    -----------
    verbose : bool
        If True, prints detailed Vina output. If False, shows only progress.
    
    Returns:
    --------
    tuple: (success, results_summary)
        - success: True if all runs completed, False if any failed
        - results_summary: Dictionary with statistics about the run
    """
    check_tools(["vina", "vina_split"])
    try:
        # Find all configuration files
        config_files = sorted(glob.glob("*_conf.txt"))
        
        if not config_files:
            print("No configuration files found in the current directory.")
            return False, {}
        
        print(f"\n{'='*70}")
        print(f"AUTODOCK VINA BATCH PROCESSING")
        print(f"{'='*70}")
        print(f"Found {len(config_files)} configuration file(s)")
        print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        # Track results
        successful = []
        failed = []
        start_time = time.time()
        
        # Process each configuration file
        for idx, config_file in enumerate(config_files, 1):
            try:
                # Get base name (removing _conf.txt)
                base_name = config_file.replace("_conf.txt", "")
                
                print(f"[{idx}/{len(config_files)}] Running AutoDock Vina: {base_name}")
                
                run_start = time.time()
                
                # Run Vina with the configuration file
                result = subprocess.run(
                    ["vina", "--config", config_file],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                run_time = time.time() - run_start
                
                # Print Vina output if verbose
                if verbose:
                    print(result.stdout)
                else:
                    # Extract just the affinity information
                    for line in result.stdout.split('\n'):
                        if 'mode' in line.lower() or '----' in line or any(char.isdigit() for char in line[:5]):
                            print(line)
                
                if result.stderr:
                    print("Warnings/Errors:")
                    print(result.stderr)
                
                print(f"✓ Completed in {run_time:.1f}s\n")
                successful.append(base_name)
                    
            except subprocess.CalledProcessError as e:
                print(f"✖ Error running Vina for {config_file}:")
                print(e.stderr)
                failed.append(config_file)
                
                # Ask if user wants to continue
                if idx < len(config_files):
                    if not prompt_yes_no("\nContinue with remaining files? (y/n): "):
                        print("Batch processing stopped by user.")
                        break
                    print()

            except Exception as e:
                print(f"✖ Unexpected error processing {config_file}:")
                print(str(e))
                failed.append(config_file)

                if idx < len(config_files):
                    if not prompt_yes_no("\nContinue with remaining files? (y/n): "):
                        print("Batch processing stopped by user.")
                        break
                    print()
        
        # Summary
        total_time = time.time() - start_time
        
        print(f"\n{'='*70}")
        print(f"BATCH PROCESSING SUMMARY")
        print(f"{'='*70}")
        print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
        print(f"Successful: {len(successful)}/{len(config_files)}")
        print(f"Failed: {len(failed)}/{len(config_files)}")
        
        if successful:
            print(f"\n✓ Successfully processed:")
            for name in successful:
                print(f"  - {name}")
        
        if failed:
            print(f"\n✖ Failed:")
            for name in failed:
                print(f"  - {name}")
        
        print(f"{'='*70}\n")
        
        results_summary = {
            'total': len(config_files),
            'successful': len(successful),
            'failed': len(failed),
            'success_rate': len(successful) / len(config_files) * 100 if config_files else 0,
            'total_time': total_time,
            'successful_files': successful,
            'failed_files': failed
        }
        
        return len(failed) == 0, results_summary
        
    except Exception as e:
        print(f"Unexpected error in run_vina: {str(e)}")
        return False, {}


def _parse_vina_affinity(pose_lines):
    """Return the affinity (kcal/mol) from a pose's REMARK VINA RESULT line, or None."""
    for line in pose_lines:
        if line.startswith("REMARK VINA RESULT:"):
            try:
                return float(line.split()[3])
            except (IndexError, ValueError):
                return None
    return None


def summarize_modes(mode_records):
    """
    Reduce per-mode docking results to the canonical-pose metrics plus a
    best-RMSD diagnostic.

    The canonical pose is ALWAYS mode 1 — Vina writes modes in descending score
    order, so mode 1 is the best-scored pose. All four downstream metrics
    (affinity, ligand RMSD, PLIF, PPS) are taken from this single pose so that
    they describe one physical binding event.

    The pose with the lowest RMSD to the reference is also identified, but only
    as a diagnostic: it lets a reviewer see when Vina's top-scored pose disagrees
    with the most reference-like pose. It is NOT used for any reported metric.
    (This is the deliberate change away from the earlier best-RMSD pose
    selection, which made the metrics describe different poses and biased the
    geometric metrics toward the reference — see METHODS_NOTES.md, limitation 1.)

    Parameters
    ----------
    mode_records : list[dict]
        One dict per mode, ordered by Vina rank, each with keys:
          "mode"     : int  (1-based Vina rank; mode 1 = best score)
          "affinity" : float | None
          "rmsd"     : float | None  (None if the pose's atoms did not match the
                                      reference, so RMSD could not be computed)

    Returns
    -------
    dict with keys: rank1_affinity, rank1_rmsd, best_rmsd_mode, best_rmsd,
        best_rmsd_mode_affinity, pose_divergence.
    """
    rank1 = next((r for r in mode_records if r["mode"] == 1), None)
    rank1_affinity = rank1["affinity"] if rank1 else None
    rank1_rmsd     = rank1["rmsd"] if rank1 else None

    best = None
    for r in mode_records:
        if r["rmsd"] is None:
            continue
        if best is None or r["rmsd"] < best["rmsd"]:
            best = r

    if best is None:
        return {
            "rank1_affinity": rank1_affinity,
            "rank1_rmsd": rank1_rmsd,
            "best_rmsd_mode": None,
            "best_rmsd": None,
            "best_rmsd_mode_affinity": None,
            "pose_divergence": False,
        }

    return {
        "rank1_affinity": rank1_affinity,
        "rank1_rmsd": rank1_rmsd,
        "best_rmsd_mode": best["mode"],
        "best_rmsd": best["rmsd"],
        "best_rmsd_mode_affinity": best["affinity"],
        "pose_divergence": best["mode"] != 1,
    }


def process_vina_output_top_mode_only(chemical, reference_pdb, ligand_resname, details_dir=None):
    """
    Keep mode 1 (Vina's best-scored pose) as the canonical pose from which ALL four
    metrics derive: the model PDB (and its downstream PLIF/PPS), the reported affinity,
    and ligand RMSD all come from this one pose, so they describe one physical state.
    vina_split produces per-mode files; per-mode RMSD vs the reference ligand is via the
    Hungarian algorithm. The lowest-RMSD-to-reference mode is recorded only as a
    diagnostic (best_rmsd_* columns + pose_divergence). Anchoring to the top-scored pose
    avoids the incoherence of best-case-per-axis metrics — see METHODS_NOTES.md.

    Returns (success, file_count, output_dir).
    """
    print(f"\n{'='*70}")
    print("PROCESSING VINA OUTPUT - RANK-1 (BEST-SCORED) POSE")
    print(f"{'='*70}\n")

    try:
        # Ensure chemical has .pdbqt extension
        if not chemical.endswith(".pdbqt"):
            chemical += ".pdbqt"

        # Remove .pdbqt for pattern matching
        ligand_id = chemical.replace(".pdbqt", "")

        # Load reference ligand coordinates (used for every protein in this batch)
        with open(reference_pdb, "r") as f:
            ref_lines = f.readlines()
        ref_coords, ref_elements = parse_hetatm_coords(ref_lines, ligand_resname)
        if len(ref_coords) == 0:
            print(f"ERROR: No heavy atoms found for ligand '{ligand_resname}' in {reference_pdb}.")
            print("Check that the ligand residue name matches the reference PDB exactly.")
            return False, 0, "docking_results"
        print(f"Reference ligand: {len(ref_coords)} heavy atoms loaded from {os.path.basename(reference_pdb)}\n")

        # Create docking_results/ as a sibling of pdbqt_files/ (one level up)
        output_folder = os.path.normpath(os.path.join(os.getcwd(), "..", "docking_results"))
        os.makedirs(output_folder, exist_ok=True)

        # Find all relevant PDBQT files
        pattern = f"*_bound_{ligand_id}.pdbqt"
        pdbqt_files = sorted(glob.glob(pattern))

        if not pdbqt_files:
            print(f"No files matching pattern '{pattern}' found.")
            print("This might mean Vina docking failed or no output was generated.")
            return False, 0, output_folder

        print(f"Found {len(pdbqt_files)} output file(s) to process")
        print(f"Using rank-1 (best-scored) pose for each file...\n")

        processed_count = 0
        score_rows = []  # accumulated per-protein data for docking_scores.csv

        for idx, pdbqt_file in enumerate(pdbqt_files, 1):
            try:
                print(f"[{idx}/{len(pdbqt_files)}] Processing: {pdbqt_file}")

                # Run vina_split to produce per-mode PDBQT files
                result = subprocess.run(
                    ["vina_split", "--input", pdbqt_file],
                    capture_output=True,
                    text=True,
                    check=True
                )
                if result.stderr:
                    print(f"  Warnings: {result.stderr.strip()}")

                base_name = os.path.splitext(pdbqt_file)[0]
                protein_name = base_name.replace(f"_bound_{ligand_id}", "")

                # ----------------------------------------------------------------
                # Per-mode scoring: for every mode whose ligand PDBQT exists collect
                #   - affinity (kcal/mol) from the REMARK VINA RESULT line
                #   - RMSD vs the reference ligand (Hungarian algorithm), or None if
                #     the pose's atom composition does not match the reference.
                # vina_split numbers modes consecutively from 1 in descending score
                # order, so the first record is always mode 1 (the best score).
                # summarize_modes() then derives the canonical (mode-1) metrics and
                # the best-RMSD diagnostic.
                # ----------------------------------------------------------------
                mode_records = []
                for mode_n in range(1, MAX_VINA_MODES + 1):
                    lig_file = f"{base_name}_ligand_{mode_n}.pdbqt"
                    if not os.path.exists(lig_file):
                        break   # stop at first gap

                    with open(lig_file, "r") as lf:
                        pose_lines = lf.readlines()

                    affinity = _parse_vina_affinity(pose_lines)
                    pose_coords, pose_elements = parse_hetatm_coords(pose_lines, ligand_resname)

                    # Validate atom composition matches reference before computing RMSD
                    pose_rmsd = None
                    if len(pose_coords) != len(ref_coords):
                        print(f"  ⚠  Mode {mode_n}: atom count mismatch "
                              f"({len(pose_coords)} vs {len(ref_coords)}) — RMSD skipped")
                    elif sorted(pose_elements) != sorted(ref_elements):
                        print(f"  ⚠  Mode {mode_n}: element mismatch — RMSD skipped")
                    else:
                        pose_rmsd = rmsd_hungarian(ref_coords, pose_coords)

                    mode_records.append({"mode": mode_n, "affinity": affinity, "rmsd": pose_rmsd})

                # Mode 1 (canonical pose) must exist; the loop starts at 1 and breaks
                # at the first gap, so an empty list means no mode 1 was produced.
                if not mode_records:
                    print(f"  ✖ No mode-1 pose found — skipping {pdbqt_file}")
                    continue

                metrics = summarize_modes(mode_records)
                rank1_affinity = metrics["rank1_affinity"]
                rank1_rmsd     = metrics["rank1_rmsd"]
                pose_divergence = metrics["pose_divergence"]

                # Report: canonical pose is mode 1; best-RMSD pose is diagnostic only
                rmsd_str = f"{rank1_rmsd:.3f} Å" if rank1_rmsd is not None else "n/a"
                print(f"  ✓ Canonical pose = mode 1 "
                      f"(affinity = {rank1_affinity} kcal/mol, RMSD to reference = {rmsd_str})")
                if pose_divergence:
                    print(f"    ⓘ Diagnostic: mode {metrics['best_rmsd_mode']} was closer to the "
                          f"reference (RMSD {metrics['best_rmsd']:.3f} Å) — recorded, NOT used")

                # Accumulate scores for docking_scores.csv. binding_affinity_kcal_mol
                # and lig_rmsd_A are the canonical (mode-1) metrics; the best_rmsd_*
                # columns are diagnostics only.
                score_rows.append({
                    "protein":                          protein_name,
                    "binding_affinity_kcal_mol":        rank1_affinity,
                    "lig_rmsd_A":                       round(rank1_rmsd, 3) if rank1_rmsd is not None else None,
                    "best_rmsd_mode":                   metrics["best_rmsd_mode"],
                    "best_rmsd_A":                      round(metrics["best_rmsd"], 3) if metrics["best_rmsd"] is not None else None,
                    "best_rmsd_mode_affinity_kcal_mol": metrics["best_rmsd_mode_affinity"],
                    "pose_divergence":                  pose_divergence,
                })

                # ----------------------------------------------------------------
                # Keep MODE 1's files (the canonical pose) with clean names. The
                # generated model PDB, and the PLIF/PPS metrics computed from it
                # downstream, all derive from this same pose.
                # ----------------------------------------------------------------
                canonical_files = {
                    f"{base_name}_ligand_1.pdbqt": f"{protein_name}_ligand.pdbqt",
                    f"{base_name}_flex_1.pdbqt":   f"{protein_name}_flex.pdbqt",
                    f"{base_name}_rigid_1.pdbqt":  f"{protein_name}_rigid.pdbqt",
                }
                moved_count = 0
                for src, dst_name in canonical_files.items():
                    if os.path.exists(src):
                        shutil.move(src, os.path.join(output_folder, dst_name))
                        moved_count += 1

                # Delete all other modes (mode-1 files are already moved out above)
                deleted_count = 0
                for i in range(1, MAX_VINA_MODES + 1):
                    for suffix in ("_ligand_", "_flex_", "_rigid_"):
                        f = f"{base_name}{suffix}{i}.pdbqt"
                        if os.path.exists(f):
                            os.remove(f)
                            deleted_count += 1

                # Move the original multi-mode bound file for record-keeping
                try:
                    shutil.move(pdbqt_file, os.path.join(output_folder, pdbqt_file))
                except OSError as e:
                    print(f"  ⚠  Could not move {pdbqt_file}: {e}")

                print(f"  ✓ Kept {moved_count} file(s) for mode 1, "
                      f"deleted {deleted_count} other mode file(s)")
                processed_count += 1

            except subprocess.CalledProcessError as e:
                print(f"  ✖ Error running vina_split: {e.stderr}")
                return False, processed_count, output_folder

            except Exception as e:
                print(f"  ✖ Unexpected error: {e}")
                return False, processed_count, output_folder

        # --------------------------------------------------------------------
        # Write docking_scores.csv so downstream tools (get_summary.py,
        # susceptibility_analysis.py) always have access to:
        #   Canonical (mode-1) metrics — all four reported metrics describe this pose:
        #     - binding_affinity_kcal_mol : mode-1 (best) docking score
        #     - lig_rmsd_A                : mode-1 pose RMSD vs reference
        #   Diagnostics (best-RMSD pose; recorded but NOT used for any metric):
        #     - best_rmsd_mode                   : which mode was closest to reference
        #     - best_rmsd_A                      : that pose's RMSD vs reference
        #     - best_rmsd_mode_affinity_kcal_mol : that pose's score
        #     - pose_divergence                  : True when best_rmsd_mode != 1
        # --------------------------------------------------------------------
        csv_dir = details_dir if details_dir and os.path.isdir(details_dir) else output_folder
        scores_csv = os.path.join(csv_dir, "docking_scores.csv")
        fieldnames = [
            "protein",
            "binding_affinity_kcal_mol",
            "lig_rmsd_A",
            "best_rmsd_mode",
            "best_rmsd_A",
            "best_rmsd_mode_affinity_kcal_mol",
            "pose_divergence",
        ]
        with open(scores_csv, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(score_rows)
        print(f"Docking scores written to: {scores_csv}")

        print(f"\n{'='*70}")
        print(f"PROCESSING COMPLETE")
        print(f"{'='*70}")
        print(f"Processed {processed_count} output file(s)")
        print(f"Results saved to: {os.path.abspath(output_folder)}")
        print(f"{'='*70}\n")

        return True, processed_count, output_folder

    except Exception as e:
        print(f"Unexpected error in process_vina_output: {str(e)}")
        return False, 0, "docking_results"


# ============================================================================
# REFERENCE PDB PARSING FUNCTIONS
# ============================================================================

def parse_reference_pdb(reference_file, ligand_resname):
    """
    Extract ligand atom ordering from the reference PDB so generated models stay
    consistent with it. Returns {'atom_order', 'atom_records', 'atom_count'}, or None.
    """
    print(f"  Parsing reference structure: {os.path.basename(reference_file)}")
    
    ligand_atoms = []
    atom_records = {}
    
    try:
        with open(reference_file, 'r') as f:
            for line in f:
                if line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    
                    if res_name == ligand_resname:
                        atom_name = line[12:16].strip()
                        ligand_atoms.append(atom_name)
                        atom_records[atom_name] = line.rstrip()
        
        if not ligand_atoms:
            print(f"  ⚠  Warning: No HETATM records found for ligand '{ligand_resname}'")
            return None
        
        print(f"  ✓ Found {len(ligand_atoms)} ligand atoms in reference")
        
        return {
            'atom_order': ligand_atoms,
            'atom_records': atom_records,
            'atom_count': len(ligand_atoms)
        }
        
    except (OSError, ValueError, IndexError) as e:
        print(f"  ✖ Error parsing reference PDB: {str(e)}")
        return None


# ============================================================================
# OBABEL-BASED MODEL GENERATION (RECOMMENDED)
# ============================================================================

def check_obabel_available():
    """
    Check if obabel is available on the system.
    
    Returns:
    --------
    bool : True if obabel is available, False otherwise
    """
    try:
        result = subprocess.run(
            ["obabel", "-V"],
            capture_output=True,
            text=True,
            check=True
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def combine_pdbqt_with_obabel(protein_name, vina_output_dir, ligand_resname, reference_info):
    """
    Combine rigid + flex + ligand PDBQT into a PDB model using obabel (rigid/flex/ligand
    converted separately, ligand atoms reordered to match the reference, then merged and
    renumbered). obabel interprets AutoDock atom types and element symbols robustly,
    which manual parsing does not. Returns the model path, or None on failure.
    """
    rigid_file = os.path.join(vina_output_dir, f"{protein_name}_rigid.pdbqt")
    flex_file = os.path.join(vina_output_dir, f"{protein_name}_flex.pdbqt")
    ligand_file = os.path.join(vina_output_dir, f"{protein_name}_ligand.pdbqt")
    
    # Check required files exist
    if not os.path.exists(rigid_file):
        print(f"  ✖ Rigid file not found: {os.path.basename(rigid_file)}")
        return None
    
    if not os.path.exists(ligand_file):
        print(f"  ✖ Ligand file not found: {os.path.basename(ligand_file)}")
        return None
    
    has_flex = os.path.exists(flex_file)
    
    # Create temporary directory for obabel conversions
    temp_dir = os.path.join(vina_output_dir, "temp_conversion")
    os.makedirs(temp_dir, exist_ok=True)
    
    try:
        # Convert ligand PDBQT to PDB
        ligand_pdb = os.path.join(temp_dir, f"{protein_name}_ligand.pdb")
        result = subprocess.run(
            ["obabel", ligand_file, "-O", ligand_pdb],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Convert rigid PDBQT to PDB
        rigid_pdb = os.path.join(temp_dir, f"{protein_name}_rigid.pdb")
        result = subprocess.run(
            ["obabel", rigid_file, "-O", rigid_pdb],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Convert flex PDBQT to PDB if present
        flex_pdb = None
        if has_flex:
            flex_pdb = os.path.join(temp_dir, f"{protein_name}_flex.pdb")
            result = subprocess.run(
                ["obabel", flex_file, "-O", flex_pdb],
                capture_output=True,
                text=True,
                check=True
            )
        
        # Read converted PDB files
        ligand_lines = []
        protein_lines = []
        flex_residues = set()
        
        # Read ligand and create ordered dictionary
        ligand_atoms = {}
        with open(ligand_pdb, 'r') as f:
            for line in f:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    atom_name = line[12:16].strip()
                    # Update residue name to match user specification
                    line = line[:17] + ligand_resname.ljust(3) + line[20:]
                    # Update residue number to a sentinel distinct from protein residues
                    line = line[:22] + LIGAND_RESNUM + line[26:]
                    ligand_atoms[atom_name] = line.rstrip()
        
        # Reorder ligand atoms to match reference structure
        if reference_info:
            for atom_name in reference_info['atom_order']:
                if atom_name in ligand_atoms:
                    # Convert ATOM to HETATM if necessary
                    line = ligand_atoms[atom_name]
                    if line.startswith("ATOM"):
                        line = "HETATM" + line[6:]
                    ligand_lines.append(line)
                else:
                    print(f"  ⚠  Warning: Atom {atom_name} from reference not found in docked ligand")
        else:
            # No reference - use order from obabel
            for atom_name, line in ligand_atoms.items():
                if line.startswith("ATOM"):
                    line = "HETATM" + line[6:]
                ligand_lines.append(line)
        
        # Read flexible residues if present
        # IMPORTANT: Track individual ATOMS, not just residues
        # Flexible files contain ONLY sidechain atoms
        # Rigid files contain backbone atoms (N, CA, C, O, H) for flexible residues
        flex_atoms = set()  # Track (res_name, res_num, atom_name) tuples
        
        if flex_pdb:
            with open(flex_pdb, 'r') as f:
                for line in f:
                    if line.startswith("ATOM"):
                        res_name = line[17:20].strip()
                        res_num = line[22:26].strip()
                        atom_name = line[12:16].strip()
                        chain = line[21] if len(line) > 21 else ' '
                        
                        # Track this specific atom as flexible
                        flex_atoms.add((chain, res_name, res_num, atom_name))
                        protein_lines.append(line.rstrip())
        
        # Read rigid protein
        # Keep ALL atoms EXCEPT those that are in the flex file
        with open(rigid_pdb, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    res_name = line[17:20].strip()
                    res_num = line[22:26].strip()
                    atom_name = line[12:16].strip()
                    chain = line[21] if len(line) > 21 else ' '
                    
                    # Only skip if this EXACT atom is in the flex file
                    if (chain, res_name, res_num, atom_name) not in flex_atoms:
                        protein_lines.append(line.rstrip())
        
        # Sort protein atoms by chain, residue number, and atom name
        def get_sort_key(line):
            chain = line[21] if len(line) > 21 else ' '
            res_num = int(line[22:26].strip()) if line[22:26].strip().lstrip('-').isdigit() else 0
            atom_name = line[12:16].strip() if len(line) > 16 else ''
            return (chain, res_num, atom_name)
        
        protein_lines.sort(key=get_sort_key)
        
        # Combine: HETATM (ligand) first, then ATOM (protein)
        all_lines = ligand_lines + protein_lines
        
        if not all_lines:
            print(f"  ✖ No valid coordinates found after conversion")
            return None
        
        # Renumber atoms sequentially
        renumbered_lines = []
        atom_counter = 1
        
        for line in all_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Replace atom number (columns 6-10, right-justified)
                new_line = line[:6] + str(atom_counter).rjust(5) + line[11:]
                renumbered_lines.append(new_line)
                atom_counter += 1
            else:
                renumbered_lines.append(line)
        
        # Add END record
        renumbered_lines.append("END")
        
        # Write output file
        output_file = os.path.join(vina_output_dir, f"{protein_name}_model.pdb")
        
        with open(output_file, 'w') as f:
            for line in renumbered_lines:
                f.write(line + '\n')
        
        # Clean up temporary directory
        shutil.rmtree(temp_dir)
        
        print(f"  ✓ Generated model with {len(ligand_lines)} ligand atoms, {len(protein_lines)} protein atoms")
        
        return output_file
        
    except subprocess.CalledProcessError as e:
        print(f"  ✖ Error running obabel: {e.stderr}")
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        return None
        
    except (OSError, ValueError, IndexError) as e:
        print(f"  ✖ Error in model generation: {str(e)}")
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
        return None


# ============================================================================
# MANUAL PDBQT-TO-PDB CONVERSION (FALLBACK)
# ============================================================================

def pdbqt_to_pdb_line(line):
    """
    Convert one PDBQT line to PDB format (fallback when obabel is unavailable):
    keeps ATOM/HETATM, strips AutoDock atom types / partial charges / activity flags /
    torsion-tree records, and drops hydrogens. Returns the line, or None to exclude it.
    Less robust than obabel and may have element-symbol issues.
    """
    # Skip hydrogen atoms
    if line.startswith("ATOM") or line.startswith("HETATM"):
        # Check element symbol (columns 77-78 in PDBQT, 76-77 in PDB)
        element = line[77:79].strip() if len(line) > 78 else line[76:78].strip()
        if element.startswith('H') or element == 'HD':
            return None
        
        # Remove activity flags (A/I) from element column
        if len(line) > 79 and line[79] in ['A', 'I']:
            line = line[:79] + ' ' + line[80:] if len(line) > 80 else line[:79]

        # Remove partial charge (columns 69-76 in PDBQT)
        if len(line) > 76:
            line = line[:69] + ' ' * 7 + line[76:]
        
        # Clean up element symbol - take first character if multi-character
        if len(element) > 1 and element[0].isalpha():
            line = line[:77] + element[0].ljust(2) + line[79:]
        
        return line
    
    # Skip AutoDock-specific records
    autodock_keywords = ["ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF", "REMARK"]
    if any(line.startswith(keyword) for keyword in autodock_keywords):
        return None
    
    # Keep other standard PDB records
    return line

def remove_hydrogens_from_pdb(pdb_file):
    """
    Strip hydrogens (by atom name or H/HD element) from a PDB file in place and
    renumber atoms — docking needs polar H, downstream analysis does not. Returns
    True on success.
    """
    try:
        # Read all lines from the PDB file
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        # Filter out hydrogen atoms
        non_hydrogen_lines = []
        
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Check atom name (columns 13-16, Python 0-indexed)
                atom_name = line[12:16].strip()
                
                # Check element symbol (columns 77-78 in PDB format, Python 0-indexed)
                element = line[76:78].strip() if len(line) > 77 else ''
                
                # Skip if atom name starts with H or element is H/HD
                if atom_name.startswith('H') or element in ['H', 'HD']:
                    continue
                
                non_hydrogen_lines.append(line)
            else:
                # Keep all non-atom records (HEADER, REMARK, END, etc.)
                non_hydrogen_lines.append(line)
        
        # Renumber atoms sequentially
        renumbered_lines = []
        atom_counter = 1
        
        for line in non_hydrogen_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Replace atom number (columns 7-11, right-justified, Python 0-indexed 6-11)
                new_line = line[:6] + str(atom_counter).rjust(5) + line[11:]
                renumbered_lines.append(new_line)
                atom_counter += 1
            else:
                renumbered_lines.append(line)
        
        # Write back to the same file
        with open(pdb_file, 'w') as f:
            for line in renumbered_lines:
                f.write(line)
        
        return True
        
    except OSError as e:
        print(f"  ✖ Error removing hydrogens from {os.path.basename(pdb_file)}: {str(e)}")
        return False

def combine_pdbqt_manual(protein_name, vina_output_dir, ligand_resname, reference_info):
    """
    Manual-parsing fallback for combine_pdbqt_with_obabel() (used when obabel is
    unavailable); may have element-symbol issues. Returns the model path, or None.
    """
    rigid_file = os.path.join(vina_output_dir, f"{protein_name}_rigid.pdbqt")
    flex_file = os.path.join(vina_output_dir, f"{protein_name}_flex.pdbqt")
    ligand_file = os.path.join(vina_output_dir, f"{protein_name}_ligand.pdbqt")
    
    if not os.path.exists(rigid_file):
        print(f"  ✖ Rigid file not found: {os.path.basename(rigid_file)}")
        return None
    
    if not os.path.exists(ligand_file):
        print(f"  ✖ Ligand file not found: {os.path.basename(ligand_file)}")
        return None
    
    has_flex = os.path.exists(flex_file)
    
    try:
        # Read ligand and create ordered dictionary
        ligand_atoms = {}
        with open(ligand_file, 'r') as f:
            for line in f:
                converted = pdbqt_to_pdb_line(line)
                if converted:
                    if converted.startswith("ATOM"):
                        converted = "HETATM" + converted[6:]
                    
                    if converted.startswith("HETATM"):
                        atom_name = converted[12:16].strip()
                        # Update residue name
                        converted = converted[:17] + ligand_resname.ljust(3) + converted[20:]
                        # Update residue number
                        converted = converted[:22] + LIGAND_RESNUM + converted[26:]
                        ligand_atoms[atom_name] = converted.rstrip()
        
        # Reorder ligand atoms to match reference
        ligand_lines = []
        if reference_info:
            for atom_name in reference_info['atom_order']:
                if atom_name in ligand_atoms:
                    ligand_lines.append(ligand_atoms[atom_name])
                else:
                    print(f"  ⚠  Warning: Atom {atom_name} from reference not found in docked ligand")
        else:
            ligand_lines = list(ligand_atoms.values())
        
        # Read flexible residues
        # IMPORTANT: Track individual ATOMS, not just residues
        flex_atoms = set()  # Track (chain, res_name, res_num, atom_name) tuples
        protein_lines = []
        
        if has_flex:
            with open(flex_file, 'r') as f:
                for line in f:
                    converted = pdbqt_to_pdb_line(line)
                    if converted and converted.startswith("ATOM"):
                        res_name = converted[17:20].strip()
                        res_num = converted[22:26].strip()
                        atom_name = converted[12:16].strip()
                        chain = converted[21] if len(converted) > 21 else ' '
                        
                        flex_atoms.add((chain, res_name, res_num, atom_name))
                        protein_lines.append(converted.rstrip())
        
        # Read rigid protein
        with open(rigid_file, 'r') as f:
            for line in f:
                converted = pdbqt_to_pdb_line(line)
                if converted and converted.startswith("ATOM"):
                    res_name = converted[17:20].strip()
                    res_num = converted[22:26].strip()
                    atom_name = converted[12:16].strip()
                    chain = converted[21] if len(converted) > 21 else ' '
                    
                    # Only skip if this EXACT atom is in flex file
                    if (chain, res_name, res_num, atom_name) not in flex_atoms:
                        protein_lines.append(converted.rstrip())
        
        # Sort protein atoms
        def get_sort_key(line):
            chain = line[21] if len(line) > 21 else ' '
            res_num = int(line[22:26].strip()) if line[22:26].strip().lstrip('-').isdigit() else 0
            atom_name = line[12:16].strip() if len(line) > 16 else ''
            return (chain, res_num, atom_name)
        
        protein_lines.sort(key=get_sort_key)
        
        # Combine and renumber
        all_lines = ligand_lines + protein_lines
        
        if not all_lines:
            print(f"  ✖ No valid coordinates found")
            return None
        
        renumbered_lines = []
        atom_counter = 1
        
        for line in all_lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                new_line = line[:6] + str(atom_counter).rjust(5) + line[11:]
                renumbered_lines.append(new_line)
                atom_counter += 1
            else:
                renumbered_lines.append(line)
        
        renumbered_lines.append("END")
        
        # Write output
        output_file = os.path.join(vina_output_dir, f"{protein_name}_model.pdb")
        
        with open(output_file, 'w') as f:
            for line in renumbered_lines:
                f.write(line + '\n')
        
        print(f"  ✓ Generated model with {len(ligand_lines)} ligand atoms, {len(protein_lines)} protein atoms")
        print(f"  ⚠  Note: Manual conversion used - element symbols may need verification")
        
        return output_file
        
    except (OSError, ValueError, IndexError) as e:
        print(f"  ✖ Error in manual model generation: {str(e)}")
        return None


# ============================================================================
# MODEL GENERATION WORKFLOW
# ============================================================================

def copy_missing_rigid_files(vina_output_dir, original_pdbqt_dir):
    """
    Copy any rigid PDBQT files still in pdbqt_files/ into docking_results/ (where the
    ligand/flex files already are). Returns the number copied.
    """
    print(f"\n  Checking for missing rigid files...")
    
    # Find all ligand files in vina_output
    ligand_files = glob.glob(os.path.join(vina_output_dir, "*_ligand.pdbqt"))
    
    if not ligand_files:
        print("    No ligand files found in docking_results")
        return 0
    
    copied_count = 0
    
    for ligand_file in ligand_files:
        # Derive the expected rigid filename
        base_name = os.path.basename(ligand_file).replace("_ligand.pdbqt", "")
        rigid_name = f"{base_name}_rigid.pdbqt"
        
        rigid_in_output = os.path.join(vina_output_dir, rigid_name)
        rigid_in_original = os.path.join(original_pdbqt_dir, rigid_name)
        
        # If rigid file missing from docking_results but exists in pdbqt_files dir
        if not os.path.exists(rigid_in_output) and os.path.exists(rigid_in_original):
            try:
                shutil.copy2(rigid_in_original, rigid_in_output)
                print(f"    ✓ Copied: {rigid_name}")
                copied_count += 1
            except OSError as e:
                print(f"    ✖ Error copying {rigid_name}: {str(e)}")
    
    if copied_count > 0:
        print(f"  Copied {copied_count} rigid file(s) to docking_results")
    else:
        print(f"  All rigid files already present (or not found in original directory)")
    
    return copied_count


def generate_models_from_pdbqt(vina_output_dir, ligand_resname, reference_file, pdbqt_dir=None):
    """
    Build a deprotonated PDB model for each docked structure: combine rigid + flex +
    ligand PDBQT (obabel if available, else manual parsing), reorder ligand atoms to the
    reference for consistent HETATM ordering, and save into models/. pdbqt_dir locates
    any missing rigid files. Returns True if at least one model was generated.
    """
    print(f"\n{'='*70}")
    print("GENERATING PDB MODELS FROM PDBQT FILES")
    print(f"{'='*70}\n")
    
    # Parse reference structure
    reference_info = parse_reference_pdb(reference_file, ligand_resname)
    
    if not reference_info:
        print("  ⚠  Warning: Could not parse reference structure")
        print("  Proceeding without reference validation")
    
    # Check obabel availability
    use_obabel = check_obabel_available()
    
    if use_obabel:
        print(f"  ✓ Open Babel detected - using obabel for conversion")
    else:
        print(f"  ⚠  Open Babel not found - using manual PDBQT parsing")
        print(f"  Note: Manual parsing may have element symbol issues")
    
    # Copy missing rigid files from pdbqt_files/ into docking_results/
    search_dir = pdbqt_dir if (pdbqt_dir and os.path.exists(pdbqt_dir)) else os.path.dirname(vina_output_dir)
    if search_dir and os.path.exists(search_dir):
        copy_missing_rigid_files(vina_output_dir, search_dir)
    
    # Find all rigid PDBQT files
    rigid_pattern = os.path.join(vina_output_dir, "*_rigid.pdbqt")
    rigid_files = sorted(glob.glob(rigid_pattern))
    
    if not rigid_files:
        print(f"\n✖ ERROR: No rigid PDBQT files found matching pattern: {rigid_pattern}")
        print(f"\nAvailable files in {vina_output_dir}:")
        all_files = os.listdir(vina_output_dir)
        for f in sorted(all_files):
            print(f"  {f}")
        return False
    
    print(f"\nFound {len(rigid_files)} structure(s) to process")
    print(f"Ligand residue name: {ligand_resname}\n")
    
    # Create models subdirectory
    models_dir = os.path.join(vina_output_dir, "models")
    os.makedirs(models_dir, exist_ok=True)
    
    successful = 0
    failed = 0
    failed_details = []
    
    # Process each structure
    for idx, rigid_file in enumerate(rigid_files, 1):
        # Extract protein name from rigid file
        protein_name = os.path.basename(rigid_file).replace("_rigid.pdbqt", "")
        
        print(f"[{idx}/{len(rigid_files)}] Processing: {protein_name}")
        
        # Generate model using appropriate method
        try:
            if use_obabel:
                model_file = combine_pdbqt_with_obabel(
                    protein_name, vina_output_dir, ligand_resname, reference_info
                )
            else:
                model_file = combine_pdbqt_manual(
                    protein_name, vina_output_dir, ligand_resname, reference_info
                )
            
            if model_file:
                # Move to models directory
                model_basename = os.path.basename(model_file)
                dest_path = os.path.join(models_dir, model_basename)
                
                shutil.move(model_file, dest_path)
                
                # Remove hydrogens from the generated model
                print(f"  Removing hydrogens from model...")
                if remove_hydrogens_from_pdb(dest_path):
                    print(f"  ✓ Model deprotonated successfully")
                else:
                    print(f"  ⚠  Warning: Could not remove hydrogens")
                
                successful += 1
            else:
                failed += 1
                failed_details.append(f"{protein_name}: Model generation returned None")
        except Exception as e:
            print(f"  ✖ Error: {str(e)}")
            failed += 1
            failed_details.append(f"{protein_name}: {str(e)}")
    
    print(f"\n{'='*70}")
    print(f"MODEL GENERATION COMPLETE")
    print(f"{'='*70}")
    print(f"Successful: {successful}/{len(rigid_files)}")
    print(f"Failed: {failed}/{len(rigid_files)}")
    
    if failed > 0 and failed_details:
        print(f"\nFailure details:")
        for detail in failed_details:
            print(f"  - {detail}")
    
    if successful > 0:
        print(f"\nModels saved to: {os.path.abspath(models_dir)}")
    print(f"{'='*70}\n")
    
    return successful > 0


# ============================================================================
# MAIN INTEGRATED WORKFLOW
# ============================================================================

def run_integrated_workflow():
    """
    End-to-end workflow: validate inputs, run Vina batch docking, process output
    (rank-1 pose), and generate reference-ordered PDB models.
    """
    print("\n" + "="*70)
    print("AUTODOCK VINA INTEGRATED WORKFLOW")
    print("Batch Docking → Top Mode Extraction → Model Generation")
    print("="*70 + "\n")

    from utils import (resolve_project_dir, load_config, save_config,
                       get_project_paths, resolve_reference_pdb)

    project_dir = resolve_project_dir()
    config = load_config(project_dir)
    paths = get_project_paths(project_dir)

    # ========================================================================
    # STEP 1: Get working directory and validate
    # ========================================================================

    pdbqt_dir = paths["pdbqt_files"]
    if not os.path.isdir(pdbqt_dir):
        print(f"ERROR: pdbqt_files/ not found at {pdbqt_dir}")
        print("Please run prep_pdbqt.py first.")
        return

    # Resolve to absolute path before chdir so we can pass it to downstream functions
    pdbqt_dir = os.path.abspath(pdbqt_dir)
    original_dir = os.getcwd()
    os.chdir(pdbqt_dir)
    print(f"Working directory: {os.getcwd()}\n")

    # ========================================================================
    # STEP 2: Validate required files for docking
    # ========================================================================

    # Check for config files
    config_count = len(glob.glob("*_conf.txt"))
    if config_count == 0:
        print("ERROR: No configuration files (*_conf.txt) found in this directory.")
        print("Please run the config file generator first.")
        os.chdir(original_dir)
        return

    print(f"Found {config_count} configuration file(s)")

    # ── Ligand name ──────────────────────────────────────────────────────────
    ligand_id = config.get("ligand_name")
    if ligand_id and os.path.exists(f"{ligand_id}.pdbqt"):
        print(f"\nUsing ligand from config: {ligand_id}.pdbqt")
    else:
        if ligand_id:
            print(f"\nWarning: {ligand_id}.pdbqt not found. Please enter manually.")
        chemical = input("\nEnter the ligand filename (with or without .pdbqt): ").strip()
        ligand_id = chemical.replace(".pdbqt", "")
        while not os.path.exists(f"{ligand_id}.pdbqt"):
            print(f"Ligand file '{ligand_id}.pdbqt' not found in directory.")
            chemical = input("Please enter the ligand filename (case sensitive): ").strip()
            ligand_id = chemical.replace(".pdbqt", "")
        config["ligand_name"] = ligand_id
        save_config(project_dir, config)

    print(f"Using ligand: {ligand_id}.pdbqt")

    # ========================================================================
    # STEP 3: Gather reference structure info and decide on workflow scope
    # ========================================================================

    print(f"\n{'='*70}")
    print("REFERENCE STRUCTURE")
    print(f"{'='*70}")
    print("\nThe reference PDB is required for two purposes:")
    print("  1. RMSD-based pose selection — the docked pose most similar to the")
    print("     reference binding mode is kept, not simply the top-scoring pose.")
    print("  2. Model generation — atom ordering and element symbols are matched")
    print("     to the reference so all output models are directly comparable.")

    # ── Ligand residue name ──────────────────────────────────────────────────
    ligand_resname = config.get("ligand_resname")
    if ligand_resname:
        print(f"\nUsing ligand residue name from config: {ligand_resname}")
    else:
        ligand_resname = input("\nEnter 3-letter ligand residue name (e.g., LIG, DHT): ").strip().upper()
        while len(ligand_resname) != 3:
            print("Ligand name must be exactly 3 characters.")
            ligand_resname = input("Enter 3-letter ligand residue name: ").strip().upper()
        config["ligand_resname"] = ligand_resname
        save_config(project_dir, config)
    print(f"Ligand residue name: {ligand_resname}")

    # ── Reference PDB path ───────────────────────────────────────────────────
    reference_file = resolve_reference_pdb(config, project_dir)
    if reference_file:
        print(f"\nUsing reference PDB from config: {os.path.basename(reference_file)}")
    else:
        reference_file = input("\nEnter path to reference PDB file: ").strip().strip('"')
        while not os.path.exists(reference_file):
            print(f"Reference file not found: {reference_file}")
            reference_file = input("Enter path to reference PDB file: ").strip().strip('"')
        config["reference_pdb"] = os.path.basename(reference_file)
        save_config(project_dir, config)
    print(f"Reference structure: {os.path.basename(reference_file)}")

    # ========================================================================
    # STEP 3b: Decide on model generation scope
    # ========================================================================

    print(f"\n{'='*70}")
    print("WORKFLOW PLANNING")
    print(f"{'='*70}")
    print("\nAfter docking completes, would you like to generate PDB models?")
    print("Model generation will:")
    print("  - Combine rigid, flex, and ligand PDBQT files")
    print("  - Convert to PDB format using Open Babel (if available)")
    print("  - Preserve ligand atom ordering from reference structure")
    print("  - Handle flexible residues automatically")
    print("  - Organize final models in a 'models' directory")

    generate_models = prompt_yes_no("\nGenerate PDB models after docking? (y/n): ")
    
    # ========================================================================
    # STEP 4: Run Vina docking
    # ========================================================================
    
    success, results = run_vina(verbose=True)
    
    if not success:
        print("\nBatch docking completed with errors.")
        print("Check the output above for details.")
        os.chdir(original_dir)
        return
    
    if results.get('successful', 0) == 0:
        print("\nNo successful docking runs to process.")
        os.chdir(original_dir)
        return
    
    # ========================================================================
    # STEP 5: Process Vina output - BEST-RMSD POSE SELECTION
    # ========================================================================

    os.makedirs(paths["results"], exist_ok=True)
    process_success, file_count, output_dir = process_vina_output_top_mode_only(
        ligand_id, reference_file, ligand_resname,
        details_dir=paths["results"]
    )
    
    if not process_success:
        print("Docking completed but there were errors processing output files.")
        os.chdir(original_dir)
        return
    
    print("Batch docking and output processing completed successfully!")
    
    # ========================================================================
    # STEP 6: Generate models if requested
    # ========================================================================
    
    if generate_models:
        model_success = generate_models_from_pdbqt(
            os.path.abspath(output_dir),
            ligand_resname,
            reference_file,
            pdbqt_dir=pdbqt_dir
        )
        
        if model_success:
            print("\n" + "="*70)
            print("WORKFLOW COMPLETED SUCCESSFULLY!")
            print("="*70)
            print(f"\nAll outputs are in: {os.path.abspath(output_dir)}")
            print(f"Final models are in: {os.path.abspath(os.path.join(output_dir, 'models'))}")
            print("="*70 + "\n")
        else:
            print("\nModel generation encountered errors. Check output above.")
    else:
        print("\n" + "="*70)
        print("WORKFLOW COMPLETED")
        print("="*70)
        print(f"\nDocking results are in: {os.path.abspath(output_dir)}")
        print("="*70 + "\n")
    
    # Modify the format of the reference file to match the generated models and copy to output directory
    if reference_file:
        hetatm_lines = []
        atom_lines = []

        with open(reference_file, "r") as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue

                # Remove hydrogens
                element = line[76:78].strip() if len(line) > 77 else ""
                if element.startswith("H"):
                    continue

                resname = line[17:20].strip()

                if line.startswith("HETATM"):
                    # Keep ONLY ligand HETATM
                    if resname != ligand_resname:
                        continue

                    # Modify residue name and number to match generated models
                    line = line[:17] + ligand_resname.ljust(3) + line[20:]
                    line = line[:22] + LIGAND_RESNUM + line[26:]

                    hetatm_lines.append(line)

                elif line.startswith("ATOM"):
                    atom_lines.append(line)

        # Combine in required order
        ordered_lines = hetatm_lines + atom_lines

        # Renumber atom serial numbers sequentially
        renumbered_lines = []
        serial = 1

        for line in ordered_lines:
            new_serial = f"{serial:5d}"
            line = line[:6] + new_serial + line[11:]
            renumbered_lines.append(line)
            serial += 1

        renumbered_lines.append("END\n")

        output_path = os.path.join(output_dir, "models", os.path.basename(reference_file))

        with open(output_path, "w") as out:
            out.writelines(renumbered_lines)

        # Check that the file was written successfully
        if os.path.exists(output_path):
            print(f"  ✓ Reference structure updated and saved to: {output_path}")
        else:
            print(f"  ✖ Error: Failed to save updated reference structure to {output_path}")

    # Return to original directory
    os.chdir(original_dir)

if __name__ == "__main__":
    run_integrated_workflow()