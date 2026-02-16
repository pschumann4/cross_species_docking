"""
Integrated PLIP Analysis and PLIF Generation Pipeline

Workflow:
1. Runs PLIP on all PDB files in specified directory
2. Processes protonated structures with OpenBabel
3. Organizes PLIP outputs (XML and protonated PDB files)
4. Automatically generates PLIFs for all test structures vs reference
5. Calculates Tanimoto similarity coefficients
6. Creates organized output with individual comparisons and summary table
"""

import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import rdkit
from rdkit import DataStructs


def euclidean3d(v1, v2):
    """
    Calculate euclidean distance for 3D coordinates.
    
    This is used for Van der Waals contact detection, where we need to identify
    protein-ligand atom pairs within a distance threshold based on their combined
    VDW radii. The explicit distance calculation allows us to filter interactions
    that PLIP might miss, particularly weak VDW contacts.
    """
    if not len(v1) == 3 and len(v2) == 3:
        return None
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def run_plip_analysis(dir_path):
    """
    Run PLIP on all PDB files in the specified directory and organize results.
    
    Parameters:
    -----------
    dir_path : str
        Path to directory containing PDB files to analyze
    
    Returns:
    --------
    str
        Path to the PLIP_results directory containing XML and protonated PDB files
    
    Methodology:
    -----------
    PLIP (Protein-Ligand Interaction Profiler) detects non-covalent interactions
    including hydrogen bonds, hydrophobic contacts, pi-stacking, salt bridges, etc.
    The -xv flag generates XML output with detailed interaction geometry, and --name
    ensures consistent file naming for downstream processing.
    
    OpenBabel is used to add hydrogens to the protonated structures, which is
    necessary for accurate representation of hydrogen bonding and protonation states
    that may vary between species or pH conditions.
    """
    print("\n" + "=" * 70)
    print("Running PLIP Analysis")
    print("=" * 70)
    
    # Change to specified directory
    os.chdir(dir_path)
    
    # Run PLIP on all PDB files
    pdb_files = [f for f in os.listdir(dir_path) if f.endswith(".pdb")]
    print(f"\nFound {len(pdb_files)} PDB files to analyze")
    
    for file in pdb_files:
        base_name = os.path.splitext(file)[0]
        try:
            print(f"  Running PLIP on {file}...")
            subprocess.run(
                ["plip", "-f", file, "-xv", "--name", base_name],
                check=True,
                capture_output=True
            )
        except subprocess.CalledProcessError as e:
            print(f"  Error processing {file}: {e}")
        except Exception as e:
            print(f"  Unexpected error processing {file}: {e}")
    
    # Process protonated files with OpenBabel
    print("\nProcessing protonated structures with OpenBabel...")
    protonated_files = [f for f in os.listdir(dir_path) if f.endswith('_protonated.pdb')]
    for file in protonated_files:
        print(f"  Adding hydrogens to {file}...")
        subprocess.run(
            ["obabel", file, "-o", "pdb", "-O", file, "-h"],
            capture_output=True
        )
    
    # Create and organize results directory
    results_dir = os.path.join(dir_path, "PLIP_results")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    print("\nOrganizing files into PLIP_results directory...")
    
    # Move XML and protonated PDB files
    files_moved = 0
    for file in os.listdir(dir_path):
        if file.endswith(".xml") or "_protonated.pdb" in file:
            src = os.path.join(dir_path, file)
            dst = os.path.join(results_dir, file)
            try:
                shutil.move(src, dst)
                files_moved += 1
            except Exception as e:
                print(f"  Warning: Could not move {file}: {e}")
    
    print(f"Moved {files_moved} files to PLIP_results folder")
    
    return results_dir


def parse_plip_xml(xml_file):
    """
    Parse PLIP XML output and return DataFrame with residue number, 
    residue type, and interaction type for each interaction.
    
    Methodology:
    -----------
    PLIP organizes interactions by type in the XML structure. We extract all
    interaction types and represent them at the residue level (not atom level)
    because residue-level patterns are more stable across species than specific
    atom positions, which can vary due to side-chain conformations.
    
    Duplicates are removed because we want to know IF a residue participates
    in each interaction type, not HOW MANY times. This creates a binary
    fingerprint that's more robust to minor structural differences that can
    exist between species.
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    data = []
    
    binding_site = root.find('bindingsite[@id="1"][@has_interactions="True"]')
    if binding_site is not None:
        interactions = binding_site.find("interactions")
        
        for interaction_type in [
            "hydrophobic_interactions",
            "hydrogen_bonds",
            "water_bridges",
            "salt_bridges",
            "pi-stacks",
            "pi_cation_interactions",
            "halogen_bonds",
            "metal_complexes",
        ]:
            interaction_elem = interactions.find(interaction_type)
            if interaction_elem is not None:
                for interaction in interaction_elem:
                    resnr = interaction.find("resnr")
                    if resnr is not None:
                        resnr = resnr.text
                        restype = interaction.find("restype").text
                        data.append({
                            "resnr": resnr,
                            "restype": restype,
                            "interaction_type": interaction_type,
                        })
    
    df = pd.DataFrame(data)
    # Remove duplicates to count interactions on a per residue basis
    df = df.drop_duplicates(keep="first")
    return df


def get_vdw_contacts(pdb_file, ligand_name):
    """
    Calculate van der Waals contacts between protein and ligand.
    
    Parameters:
    -----------
    pdb_file : str
        Path to protonated PDB file
    ligand_name : str
        Ligand identifier as found in HETATM records
    
    Returns:
    --------
    DataFrame
        Van der Waals contacts with residue number, type, and interaction type
    
    Methodology:
    -----------
    VDW contacts are detected by comparing inter-atomic distances to the sum of
    VDW radii plus a 0.6 Å tolerance. This captures weak contacts that may not
    form classical interaction types but still contribute to binding affinity.
    
    This is particularly important for cross-species comparisons because even if
    specific H-bonds or salt bridges differ due to sequence variation, conserved
    VDW contacts can indicate similar binding modes. The tolerance of 0.6 Å
    accounts for thermal fluctuations and coordinate uncertainty.
    """
    vdw_radii = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}
    protein_coordinates = []
    ligand_coordinates = []
    
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom = line[12:16].strip()[0]
                resnr = line[22:26].strip()
                restype = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                protein_coordinates.append([atom, resnr, restype, x, y, z])
            elif line.startswith("HETATM") and ligand_name in line:
                atom = line[12:16].strip()[0]
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ligand_coordinates.append([atom, x, y, z])
    
    if protein_coordinates and ligand_coordinates:
        all_rows = []
        for lig_coord in ligand_coordinates:
            for prot_coord in protein_coordinates:
                dist = euclidean3d(
                    [lig_coord[1], lig_coord[2], lig_coord[3]],
                    [prot_coord[3], prot_coord[4], prot_coord[5]]
                )
                all_rows.append({
                    "HETATM": lig_coord[0],
                    "LIG x": lig_coord[1],
                    "LIG y": lig_coord[2],
                    "LIG z": lig_coord[3],
                    "ATOM": prot_coord[0],
                    "resnr": prot_coord[1],
                    "restype": prot_coord[2],
                    "PROT x": prot_coord[3],
                    "PROT y": prot_coord[4],
                    "PROT z": prot_coord[5],
                    "DIST": dist
                })
        
        df = pd.DataFrame(all_rows)
        
        # Calculate vdw_radii and interactions
        df["vdw_radii"] = df.apply(
            lambda row: vdw_radii.get(row["HETATM"][0], 0) + 
                       vdw_radii.get(row["ATOM"][0], 0) + 0.6,
            axis=1
        )
        df["vdw_interaction"] = df["DIST"] < df["vdw_radii"]
        
        vdw_df = df[df["vdw_interaction"]].copy()
        vdw_df = vdw_df.groupby(["resnr", "restype"]).first().reset_index()
        vdw_df = vdw_df[["resnr", "restype"]].assign(interaction_type="vdw_contact")
    else:
        vdw_df = pd.DataFrame(columns=["resnr", "restype", "interaction_type"])
    
    return vdw_df


def merge_plifs(ref_df, test_df):
    """
    Merge reference and test PLIFs to create comparison fingerprint.
    
    Parameters:
    -----------
    ref_df : DataFrame
        Reference PLIF with resnr, restype, interaction_type
    test_df : DataFrame
        Test PLIF with resnr, restype, interaction_type
    
    Returns:
    --------
    DataFrame
        Merged fingerprint with binary columns indicating presence in ref/test
    
    Methodology:
    -----------
    The outer merge creates a union of all interactions seen in either structure,
    allowing us to identify both shared and unique interactions. This is critical
    for cross-species analysis because we need to know:
    1. Which interactions are conserved (present in both)
    2. Which are lost in the test species (only in reference)
    3. Which are gained in the test species (only in test)
    
    The binary encoding (0/1 for absent/present) enables Tanimoto similarity
    calculation, which properly accounts for both presence and absence of features.
    """
    ref_copy = ref_df.copy().rename(columns={
        "resnr": "resnr_ref",
        "restype": "restype_ref",
        "interaction_type": "interaction_type_ref"
    })
    
    test_copy = test_df.copy().rename(columns={
        "resnr": "resnr_test",
        "restype": "restype_test",
        "interaction_type": "interaction_type_test"
    })
    
    # Merge on matching interactions
    merged_df = pd.merge(
        ref_copy,
        test_copy,
        left_on=["resnr_ref", "restype_ref", "interaction_type_ref"],
        right_on=["resnr_test", "restype_test", "interaction_type_test"],
        how="outer"
    )
    
    # Create binary columns for presence in ref and test
    merged_df = merged_df.assign(
        ref=merged_df["resnr_ref"].notnull().astype(int),
        test=merged_df["resnr_test"].notnull().astype(int)
    )
    
    # Fill missing values and clean up
    merged_df["resnr_ref"] = merged_df["resnr_ref"].fillna(merged_df["resnr_test"])
    merged_df["restype_ref"] = merged_df["restype_ref"].fillna(merged_df["restype_test"])
    merged_df["interaction_type_ref"] = merged_df["interaction_type_ref"].fillna(
        merged_df["interaction_type_test"]
    )
    
    # Final DataFrame with clean column names
    result_df = merged_df[["resnr_ref", "restype_ref", "interaction_type_ref", "ref", "test"]].rename(
        columns={
            "resnr_ref": "resnr",
            "restype_ref": "restype",
            "interaction_type_ref": "interaction"
        }
    )
    
    return result_df


def calculate_tanimoto(plif_df):
    """
    Calculate Tanimoto similarity coefficient from merged PLIF DataFrame.
    
    Methodology:
    -----------
    Tanimoto coefficient (also called Jaccard index for binary data) measures
    similarity as: (shared features) / (total unique features)
    
    The coefficient ranges from 0 (no overlap) to 1 (perfect match), providing
    an intuitive measure of binding site similarity across species.
    """
    # Convert binary columns to bit strings
    bit_string_ref = str(plif_df["ref"].to_numpy())
    bit_string_test = str(plif_df["test"].to_numpy())
    
    # Create RDKit bit vectors
    plif_ref = rdkit.DataStructs.cDataStructs.CreateFromBitString(bit_string_ref)
    plif_test = rdkit.DataStructs.cDataStructs.CreateFromBitString(bit_string_test)
    
    # Calculate Tanimoto coefficient
    tanimoto = DataStructs.TanimotoSimilarity(plif_ref, plif_test)
    return round(tanimoto, 3)


def generate_plifs(plip_results_dir, ligand_name, ref_pdb=None):
    """
    Generate PLIFs for all test structures compared to reference.
    
    Parameters:
    -----------
    plip_results_dir : str
        Path to directory containing PLIP XML and protonated PDB files
    ligand_name : str
        Ligand identifier as found in HETATM records
    ref_pdb : str, optional
        Name of reference PDB file (if None, will prompt user)
    
    Returns:
    --------
    DataFrame
        Summary table with Tanimoto similarities for all comparisons
    
    Methodology:
    -----------
    This function implements the core PLIF comparison workflow.
    
    Individual PLIF comparison files are saved for further inspection.
    """
    print("\n" + "=" * 70)
    print("Generating Protein-Ligand Interaction Fingerprints")
    print("=" * 70)
    
    os.chdir(plip_results_dir)
    
    # Find reference PDB file if not specified
    if ref_pdb is None:
        pdb_files = [i for i in os.listdir(plip_results_dir) 
                     if i.endswith("_protonated.pdb")]
        
        # Try to auto-detect reference file
        ref_candidates = [f for f in pdb_files if f.startswith("ref_")]
        
        if len(ref_candidates) == 1:
            ref = input(f"\nIs {ref_candidates[0]} the reference PDB file? (y/n): ")
            ref = ref.lower()
            while ref not in ["y", "n"]:
                ref = input("Please enter y or n: ")
            if ref == "y":
                ref_pdb = ref_candidates[0]
        
        if ref_pdb is None:
            print("\nAvailable protonated PDB files:")
            for i, file in enumerate(pdb_files, 1):
                print(f"  {i}. {file}")
            ref_pdb = input("\nEnter the name of the protonated reference PDB file: ")
            if not ref_pdb.endswith(".pdb"):
                ref_pdb += ".pdb"
            while not os.path.exists(ref_pdb):
                print('Error: This file does not appear to exist.')
                ref_pdb = input("Enter the file name of the protonated reference PDB file: ")
                if not ref_pdb.endswith(".pdb"):
                    ref_pdb += ".pdb"
    
    # Get reference base name
    ref_pdb_name = ref_pdb.replace("_protonated.pdb", "")
    ref_xml = ref_pdb_name + ".xml"
    
    print(f"\nProcessing reference structure: {ref_pdb_name}")
    print("=" * 70)
    
    # Parse reference PLIP data
    print(f"Parsing PLIP XML report for {ref_pdb_name}...")
    ref_plip_df = parse_plip_xml(ref_xml)
    
    # Get reference Van der Waals contacts
    print(f"Calculating Van der Waals contacts for {ref_pdb_name}...")
    ref_vdw_df = get_vdw_contacts(ref_pdb, ligand_name)
    
    # Combine reference PLIP and VDW data
    ref_df = pd.concat([ref_plip_df, ref_vdw_df]).reset_index(drop=True)
    
    print(f"\nReference PLIF contains {len(ref_df)} interactions")
    print(f"  - PLIP-detected: {len(ref_plip_df)}")
    print(f"  - VDW contacts: {len(ref_vdw_df)}")
    
    print("\nProcessing test structures...")
    print("=" * 70)
    
    # Process all test structures
    results = []
    
    for pdb in os.listdir(plip_results_dir):
        if pdb.endswith("_protonated.pdb") and pdb != ref_pdb:
            pdb_name = pdb.replace("_protonated.pdb", "")
            xml = pdb_name + ".xml"
            
            print(f"\nProcessing {pdb_name}...")
            
            # Parse test PLIP data
            test_plip_df = parse_plip_xml(xml)
            
            # Get test Van der Waals contacts
            test_vdw_df = get_vdw_contacts(pdb, ligand_name)
            
            # Combine test PLIP and VDW data
            test_df = pd.concat([test_plip_df, test_vdw_df]).reset_index(drop=True)
            
            print(f"  Test PLIF contains {len(test_df)} interactions")
            print(f"    - PLIP-detected: {len(test_plip_df)}")
            print(f"    - VDW contacts: {len(test_vdw_df)}")
            
            # Merge reference and test PLIFs
            merged_plif = merge_plifs(ref_df, test_df)
            
            # Calculate Tanimoto coefficient
            tanimoto = calculate_tanimoto(merged_plif)
            
            shared = merged_plif[(merged_plif["ref"] == 1) & (merged_plif["test"] == 1)].shape[0]
            print(f"  Shared interactions: {shared}")
            print(f"  Tanimoto similarity: {tanimoto}")
            
            # Save individual PLIF comparison
            output_file = f"{pdb_name}_vs_{ref_pdb_name}_PLIF.txt"
            with open(output_file, "w") as f:
                f.write(merged_plif.to_string())
            
            # Add to results
            results.append({
                "test_structure": pdb_name,
                "reference_structure": ref_pdb_name,
                "tanimoto_similarity": tanimoto,
                "ref_interactions": len(ref_df),
                "test_interactions": len(test_df),
                "shared_interactions": shared
            })
    
    # Create results summary table
    print("\n" + "=" * 70)
    print("Creating summary table...")
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("tanimoto_similarity", ascending=False)
    
    # Save summary table
    summary_file = "plif_similarity_summary.csv"
    results_df.to_csv(summary_file, index=False)
    
    print(f"\nSummary table saved to: {summary_file}")
    
    # Organize output files
    print("\n" + "=" * 70)
    print("Organizing output files...")
    print("=" * 70)
    
    output_dir = os.path.join(plip_results_dir, "PLIF_files")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Move PLIF files to output directory
    files_moved = 0
    for file in os.listdir(plip_results_dir):
        if file.endswith("_PLIF.txt") or file == summary_file:
            try:
                shutil.move(file, os.path.join(output_dir, file))
                files_moved += 1
            except shutil.Error:
                print(f"  Warning: {file} already exists in output directory")
    
    print(f"Moved {files_moved} files to PLIF_files folder")
    
    return results_df


def main():
    """
    Main integrated workflow for PLIP analysis and PLIF generation.
    
    Workflow:
    1. User provides directory with PDB files (reference + test structures)
    2. PLIP analysis runs on all structures, generating XML reports
    3. Protonated structures are processed with OpenBabel
    4. Results are organized into plip_results directory
    5. User specifies ligand ID and confirms reference structure
    6. PLIFs are generated for all test structures vs reference
    7. Summary table with Tanimoto similarities is created
    8. Individual comparison files are organized for detailed inspection
    
    Benefits of integration:
    - Reduces manual steps and potential for user error
    - Maintains consistent directory structure throughout workflow
    - Provides continuous progress feedback across both stages
    - Ensures all intermediate files are properly organized
    - Allows for error handling across the entire pipeline
    """
    print("\n" + "=" * 70)
    print("PLIP ANALYSIS AND PLIF GENERATION PIPELINE")
    print("=" * 70)
    print("\nThis script will:")
    print("1. Run PLIP analysis on all PDB models in your directory")
    print("2. Generate protein-ligand interaction fingerprints (PLIFs)")
    print("3. Compare all test structures against the reference structure")
    print("4. Calculate Tanimoto similarity coefficients")
    print("5. Create organized output files for analysis")
    
    # Get directory path from user
    dir_path = input("\nEnter the path to the folder containing docked PDB models: ")
    while not os.path.exists(dir_path):
        print("Error: Directory does not exist.")
        dir_path = input("Please enter a valid directory path: ")
    
    # Run PLIP analysis
    plip_results_dir = run_plip_analysis(dir_path)
    
    # Get ligand name from user
    ligand = input("\nEnter the ligand ID as it appears in the PDB models (e.g., 'UNL', 'LIG'): ")
    
    # Generate PLIFs
    results_df = generate_plifs(plip_results_dir, ligand)
    
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"\nResults are organized in: {plip_results_dir}")
    print(f"  - XML reports and protonated PDBs in main directory")
    print(f"  - PLIF comparisons and summary in PLIF_files subdirectory")


if __name__ == "__main__":
    main()