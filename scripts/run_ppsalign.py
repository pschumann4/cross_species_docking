"""
Integrated workflow for binding pocket extraction and structural alignment.

The functions in this script perform the following steps:
1. Extract binding pocket structures from protein-ligand complexes
2. Perform structural alignment using PPSalign to compare pockets

Workflow:
    - Extract binding pockets from a directory of PDB files
    - Identify a reference pocket structure
    - Compare all extracted pockets to the reference using PPSalign
    - Output similarity scores for cross-species comparison

"""

import os
import sys
import shutil
import subprocess
import numpy as np
from pathlib import Path


# ============================================================================
# BINDING POCKET EXTRACTION FUNCTIONS (from get_pocket.py)
# ============================================================================

def euclidean3d(v1, v2):
    """
    Faster implementation of euclidean distance for the 3D case.
    """
    if not len(v1) == 3 and len(v2) == 3:
        return None
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def centroid(coo):
    """
    Calculates the centroid from a 3D point cloud and returns the coordinates.
    
    Parameters:
        coo: Array of coordinate arrays
    
    Returns:
        centroid coordinates as list
    """
    return list(
        map(
            np.mean,
            (([c[0] for c in coo]), ([c[1] for c in coo]), ([c[2] for c in coo])),
        )
    )


def min_dist(pdb_file, ligand):
    """
    Determines the residues within 7.5 Angstroms of the ligand in the PDB file.
    
    This uses an all-atom distance calculation approach where the minimum
    distance between any ligand atom and any residue atom is computed.
    
    Parameters:
        pdb_file: Path to PDB file
        ligand: Three-letter ligand residue name
    
    Returns:
        Dictionary mapping residue numbers to minimum distances
    """
    ligand_coords = []
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    
    for line in lines:
        if line.startswith("HETATM"):
            if line[17:20].strip() == ligand:
                ligand_coords.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
    
    res_coords = {}
    for line in lines:
        if line.startswith("ATOM"):
            resnr = line[22:26]
            coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if resnr not in res_coords:
                res_coords[resnr] = [coords]
            else:
                res_coords[resnr].append(coords)
    
    min_dist_dict = {}
    for resnr, coords in res_coords.items():
        min_dist_dict[resnr] = 100
        for ligand_coord in ligand_coords:
            for coord in coords:
                dist = euclidean3d(ligand_coord, coord)
                if dist < min_dist_dict[resnr]:
                    min_dist_dict[resnr] = dist
    
    min_dist_dict = {resnr: dist for resnr, dist in min_dist_dict.items() if dist < 7.5}
    return min_dist_dict


def get_bindingsite(pdb_file, ligand, verbose=True):
    """
    Determines which residues in the PDB model are within 7.5 Angstroms of the ligand.
    
    Uses a two-step approach:
    1. Centroid-based distance calculation with ligand radius correction
    2. All-atom minimum distance verification
    
    Parameters:
        pdb_file: Path to PDB file
        ligand: Three-letter ligand residue name
        verbose: Whether to print progress information
    
    Returns:
        List of residue numbers in the binding site (sorted)
    """
    ligand_coords = []
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    
    for line in lines:
        if line.startswith("HETATM"):
            if line[17:20].strip() == ligand:
                ligand_coords.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
    
    if not ligand_coords:
        if verbose:
            print(f"Ligand {ligand} not found in PDB file {pdb_file}.")
        return None
    
    ligand_centroid = centroid(ligand_coords)
    if verbose:
        print(f"Ligand centroid: {ligand_centroid}")
    
    res_coords = {}
    for line in lines:
        if line.startswith("ATOM"):
            resnr = line[22:26]
            coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            if resnr not in res_coords:
                res_coords[resnr] = [coords]
            else:
                res_coords[resnr].append(coords)
    
    res_centroids = {}
    for resnr, coords in res_coords.items():
        res_centroids[resnr] = centroid(coords)
    
    BS_DIST = 7.5
    max_dist = 0
    for coords in ligand_coords:
        dist = euclidean3d(ligand_centroid, coords)
        if dist > max_dist:
            max_dist = dist
    
    cutoff = BS_DIST + max_dist
    if verbose:
        print(f"Cutoff distance: {cutoff}")
    
    bindingsite_resnr = []
    for resnr, coords in res_centroids.items():
        if euclidean3d(ligand_centroid, coords) < cutoff:
            bindingsite_resnr.append(resnr)
    
    bindingsite_resnr = list(set(bindingsite_resnr))
    min_dist_res = min_dist(pdb_file, ligand)
    bindingsite_resnr = [i for i in bindingsite_resnr if i in min_dist_res]
    bindingsite_resnr.sort()
    
    if verbose:
        print(f"Number of residues in the binding site: {len(bindingsite_resnr)}")
    
    return bindingsite_resnr


def get_pdb_bindingsite(pdb_file, bindingsite_res):
    """
    Extracts the binding site residues from the PDB file and creates two output files:
    1. A PDB file containing only the binding site residues
    2. A .poc file (pocket coordinate file) for PPSalign
    
    The .poc format requires:
    - Header line: POC [pocket_name]
    - ATOM records for all atoms in binding site residues (PDB format)
    - TER line at the end
    
    Parameters:
        pdb_file: Path to input PDB file
        bindingsite_res: List of residue numbers in the binding site
    
    Returns:
        Tuple of (bindingsite_pdb_path, poc_file_path)
    """
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    
    # Create output file paths
    base_name = os.path.basename(pdb_file).split(".pdb")[0]
    output_dir = os.path.dirname(pdb_file)
    bindingsite_pdb = os.path.join(output_dir, f"{base_name}_bindingsite.pdb")
    poc_file = os.path.join(output_dir, f"{base_name}_bindingsite.poc")
    
    # Extract binding site residues
    bindingsite_lines = []
    for line in lines:
        if line.startswith("ATOM"):
            resnr = line[22:26]
            if resnr in bindingsite_res:
                bindingsite_lines.append(line)
    
    # Write PDB file with TER record
    with open(bindingsite_pdb, "w") as f:
        for line in bindingsite_lines:
            f.write(line)
        f.write("TER\n")
    
    # Write POC file with required format:
    # POC [pocket_name]
    # [ATOM records]
    # TER
    with open(poc_file, "w") as f:
        # Write POC header with pocket name
        f.write(f"POC {base_name}\n")
        # Write all ATOM records from binding site
        for line in bindingsite_lines:
            f.write(line)
        # Write TER line
        f.write("TER\n")
    
    return bindingsite_pdb, poc_file


def extract_pockets_from_directory(pdb_dir, ligand, output_subdir="binding_sites", verbose=True):
    """
    Extract binding pockets from all PDB files in a directory.
    
    Creates .pdb and .poc files for each binding pocket.
    The .poc format follows PPSalign requirements:
    - POC [name] header
    - All atoms from binding site residues in PDB format
    - TER terminator
    """
    # Create output directory
    output_dir = os.path.join(pdb_dir, output_subdir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        if verbose:
            print(f"Created output directory: {output_dir}")
    
    # Get list of PDB files
    pdb_files = [
        os.path.join(pdb_dir, f) for f in os.listdir(pdb_dir) 
        if f.endswith(".pdb") and "bindingsite" not in f
    ]
    
    if not pdb_files:
        print(f"No PDB files found in {pdb_dir}")
        return {}
    
    pocket_files = {}
    
    for pdb_file in pdb_files:
        base_name = os.path.basename(pdb_file).split(".pdb")[0]
        
        if verbose:
            print(f"\nProcessing {base_name}...")
        
        # Read file to check for ligand and protonation
        with open(pdb_file, "r") as f:
            lines = f.readlines()
        
        # Check if ligand is present
        if ligand not in [line[17:20].strip() for line in lines if line.startswith("HETATM")]:
            if verbose:
                print(f"  Ligand {ligand} not found. Skipping.")
            continue
        
        # Check if file is protonated (contains hydrogen atoms)
        if any(line[12:16].strip().startswith("H") for line in lines if line.startswith("ATOM")):
            if verbose:
                print(f"  File appears to be protonated. Skipping.")
            continue
        
        # Extract binding site
        bindingsite_res = get_bindingsite(pdb_file, ligand, verbose=verbose)
        
        if bindingsite_res is None:
            continue
        
        # Define output file paths
        temp_pdb = os.path.join(output_dir, f"{base_name}_bindingsite.pdb")
        temp_poc = os.path.join(output_dir, f"{base_name}_bindingsite.poc")
        
        # Extract binding site residues (all atoms from those residues)
        bindingsite_lines = []
        
        for line in lines:
            if line.startswith("ATOM"):
                resnr = line[22:26]
                if resnr in bindingsite_res:
                    bindingsite_lines.append(line)
        
        # Write PDB file
        with open(temp_pdb, "w") as f:
            f.writelines(bindingsite_lines)
            f.write("TER\n")
        
        # Write POC file in PPSalign format
        with open(temp_poc, "w") as f:
            # POC header with pocket name
            f.write(f"POC {base_name}\n")
            # All ATOM records from binding site
            f.writelines(bindingsite_lines)
            # TER terminator
            f.write("TER\n")
        
        # Verify both files exist
        if os.path.exists(temp_pdb) and os.path.exists(temp_poc):
            pocket_files[base_name] = (temp_pdb, temp_poc)
            if verbose:
                print(f"  Created: {os.path.basename(temp_pdb)}")
                print(f"  Created: {os.path.basename(temp_poc)} (POC format)")
        else:
            if verbose:
                print(f"  ERROR: Failed to create output files!")
    
    if verbose:
        print(f"\nExtracted {len(pocket_files)} binding pockets to {output_dir}")
    
    return pocket_files


# ============================================================================
# PPSALIGN STRUCTURAL COMPARISON FUNCTIONS (from run_ppsalign.py)
# ============================================================================

def run_ppsalign(query_poc, template_poc, output_file):
    """
    Run PPSalign to compare two binding pocket structures.
    
    PPSalign (Protein Pocket Structure Alignment) compares binding sites
    based on their 3D structure, providing a similarity score.
    
    Parameters:
        query_poc: Path to query .poc file
        template_poc: Path to template/reference .poc file
        output_file: Path to save PPSalign output
    
    Returns:
        True if successful, False otherwise
    """
    try:
        with open(output_file, "w") as f:
            subprocess.run(
                ["PPSalign", query_poc, template_poc],
                stdout=f,
                stderr=subprocess.PIPE
            )
        return True
    except FileNotFoundError:
        print("Error: PPSalign command not found. Make sure it's in your PATH.")
        return False
    except Exception as e:
        print(f"Error running PPSalign: {e}")
        return False


def compare_pockets_to_reference(poc_dir, template_poc, output_subdir="PPS_files", verbose=True):
    """
    Compare all binding pockets in a directory to a reference pocket.
    
    This function runs PPSalign for each pocket structure, comparing them
    to a validated reference structure (typically from a species with known
    susceptibility).
    
    Parameters:
        poc_dir: Directory containing .poc files
        template_poc: Name of the reference .poc file (should be in poc_dir)
        output_subdir: Name of subdirectory to create for output files
        verbose: Whether to print progress information
    
    Returns:
        Dictionary mapping base filenames to output file paths
    """
    # Create output directory
    output_dir = os.path.join(poc_dir, output_subdir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        if verbose:
            print(f"Created output directory: {output_dir}")
    
    # Get template path
    template_path = os.path.join(poc_dir, template_poc)
    if not os.path.exists(template_path):
        print(f"Error: Template file {template_poc} not found in {poc_dir}")
        return {}
    
    # Get all .poc files
    poc_files = [f for f in os.listdir(poc_dir) if f.endswith(".poc")]
    
    if not poc_files:
        print(f"No .poc files found in {poc_dir}")
        return {}
    
    if verbose:
        print(f"\nComparing {len(poc_files) - 1} pockets to reference: {template_poc}")
    
    results = {}
    
    for poc_file in poc_files:
        # Skip the template file itself
        if poc_file == template_poc:
            continue
        
        base_name = poc_file.split(".poc")[0]
        query_path = os.path.join(poc_dir, poc_file)
        output_file = os.path.join(output_dir, f"{base_name}_PPS.txt")
        
        if verbose:
            print(f"  Calculating PPS-Score for {poc_file}...")
        
        success = run_ppsalign(query_path, template_path, output_file)
        
        if success:
            results[base_name] = output_file
    
    if verbose:
        print(f"\nPPS-Score calculation completed. Results saved to {output_dir}")
    
    return results


# ============================================================================
# INTEGRATED WORKFLOW FUNCTIONS
# ============================================================================

def identify_reference_pocket(poc_files, auto_detect=True):
    """
    Identify the reference/template pocket structure.
    
    Parameters:
        poc_files: List of .poc filenames
        auto_detect: Whether to automatically detect files starting with "ref_"
    
    Returns:
        Name of the reference .poc file
    """
    # Auto-detect reference files
    if auto_detect:
        ref_files = [f for f in poc_files if f.startswith("ref_")]
        if len(ref_files) == 1:
            print(f"Auto-detected reference: {ref_files[0]}")
            return ref_files[0]
        elif len(ref_files) > 1:
            print(f"Multiple potential reference files found: {ref_files}")
    
    # Manual selection
    print("\nAvailable .poc files:")
    for i, poc in enumerate(poc_files, 1):
        print(f"  {i}. {poc}")
    
    while True:
        choice = input("\nEnter the number of the reference/template file: ")
        try:
            idx = int(choice) - 1
            if 0 <= idx < len(poc_files):
                return poc_files[idx]
            else:
                print("Invalid number. Please try again.")
        except ValueError:
            print("Please enter a number.")


def run_integrated_workflow(pdb_dir, ligand, verbose=True):
    """
    Complete integrated workflow: pocket extraction → structural alignment.
    
    This workflow:
    1. Extracts binding pockets from all PDB files in the directory
    2. Identifies or prompts for a reference pocket structure
    3. Compares all pockets to the reference using PPSalign
    4. Outputs similarity scores for cross-species comparison
    
    Parameters:
        pdb_dir: Directory containing PDB files
        ligand: Three-letter ligand residue name
        verbose: Whether to print progress information
    
    Returns:
        Dictionary containing:
            - 'pocket_files': Mapping of filenames to (pdb, poc) paths
            - 'pps_results': Mapping of filenames to PPS result paths
            - 'reference': Name of the reference structure used
    """
    
    # Step 1: Extract binding pockets
    if verbose:
        print("\n" + "=" * 70)
        print("EXTRACTING BINDING POCKETS")
        print("=" * 70)
    
    pocket_files = extract_pockets_from_directory(pdb_dir, ligand, verbose=verbose)
    
    if not pocket_files:
        print("No pockets extracted. Exiting.")
        return None
    
    # Step 2: Identify reference pocket
    if verbose:
        print("\n" + "=" * 70)
        print("IDENTIFYING REFERENCE POCKET")
        print("=" * 70)
    
    poc_dir = os.path.join(pdb_dir, "binding_sites")
    poc_files = [f for f in os.listdir(poc_dir) if f.endswith(".poc")]
    
    # Identify reference pocket (auto-detect)
    template_poc = identify_reference_pocket(poc_files, auto_detect=True)
    
    # Step 3: Run PPSalign comparisons
    if verbose:
        print("\n" + "=" * 70)
        print("RUNNING PPSALIGN COMPARISONS")
        print("=" * 70)
    
    pps_results = compare_pockets_to_reference(poc_dir, template_poc, verbose=verbose)
    
    # Summary
    if verbose:
        print("\n" + "=" * 70)
        print("WORKFLOW COMPLETE")
        print("=" * 70)
        print(f"\nExtracted pockets: {len(pocket_files)}")
        print(f"Reference structure: {template_poc}")
        print(f"PPS comparisons completed: {len(pps_results)}")
        print(f"\nResults location: {os.path.join(poc_dir, 'PPS_files')}")
    
    return {
        'pocket_files': pocket_files,
        'pps_results': pps_results,
        'reference': template_poc
    }


# ============================================================================
# COMMAND LINE INTERFACE
# ============================================================================

def main():
    """
    Command line interface for the integrated workflow.
    """
    print("\n" + "=" * 70)
    print("BINDING POCKET EXTRACTION AND STRUCTURAL ALIGNMENT")
    print("=" * 70)
    
    # Get input directory
    while True:
        pdb_dir = input("Enter the directory containing PDB files: ").strip()
        if os.path.exists(pdb_dir) and os.path.isdir(pdb_dir):
            break
        print(f"Directory '{pdb_dir}' not found. Please try again.")
    
    # Get ligand name
    ligand = input("Enter the ligand ID (3-letter code from PDB): ").strip()
    
    # Run workflow
    results = run_integrated_workflow(
        pdb_dir=pdb_dir,
        ligand=ligand,
        verbose=True
    )


if __name__ == "__main__":
    main()