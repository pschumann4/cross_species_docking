"""
 This script is used for generating Protein-Ligand Interaction Fingerprints (PLIFs) 
 from the .xml output of the Protein-Ligand Interaction Profiler (PLIP) tool.                              
 The PLIFs include all intermolecular interactions counted by PLIP as well as 
 van der Waals contacts. This script will also calculate the Tanimoto coefficients 
 for the inputted PDBs. Each PLIF will be saved as a .txt file along with a heatmap 
 of the PLIF matrix. To execute the script, type 'python plif_matrix.py' at the 
 command line and follow the prompts.  
"""

import os
import shutil
import time
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdkit
import seaborn as sns
from rdkit import DataStructs
from rdkit.DataStructs import TanimotoSimilarity


def euclidean3d(v1, v2):
    """
    Faster implementation of euclidean distance for the 3D case.
    """
    if not len(v1) == 3 and len(v2) == 3:
        return None
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def get_bindingsite_info(xml_file):
    """
    Parses the PLIP XML output file and returns a data frame with the residue
    number, residue type, and whether or not the residue is in contact with the
    ligand. These residues are considered as those that define the binding site.
    """
    # Parse the XML file
    tree = ET.parse(xml_file)
    # Get the root element
    root = tree.getroot()
    # Find the binding site element with the specified attribute values
    binding_site = root.find('bindingsite[@id="1"][@has_interactions="True"]')
    # Make sure the binding site element is not None
    if binding_site is not None:
        # Find the bs_residues element
        bs_residues = binding_site.find("bs_residues")
        # Find all bs_residue elements
        residues = bs_residues.findall("bs_residue")
        # Create an empty df to store residue info
        data = []
        # Get residue names
        for residue in residues:
            # Get the 'aa' attribute
            aa = residue.attrib["aa"]
            # Get residue names
            res_id = residue.text
            # Get residue contacts
            contact = residue.attrib["contact"]
            data.append({"resnr": res_id, "restype": aa, "contact": contact})
    # Create data frame
    df = pd.DataFrame(data)
    # Remove chain ID from res ID column
    df["resnr"] = df["resnr"].str.replace(r"\D", "", regex=True)
    # Return the df
    return df


def parse_plip_xml(xml_file):
    """
    Parses the PLIP XML output file and returns a data frame with the
    residue number, residue type, and interaction type for each interaction.
    This is used to generate the PLIF for all interactions other than
    van der Waals contacts.
    """
    # Parse the XML file into an ElementTree object
    tree = ET.parse(xml_file)
    # Get the root element of the tree
    root = tree.getroot()
    # Initialize an empty list to store the data
    data = []
    # Get the binding site element
    binding_site = root.find('bindingsite[@id="1"][@has_interactions="True"]')
    # Make sure the binding site element is not None
    if binding_site is not None:
        # Get the interactions element
        interactions = binding_site.find("interactions")
        # Iterate over the specified interaction types
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
            # Get the interactions element for the current interaction type
            interaction_elem = interactions.find(interaction_type)
            if interaction_elem is not None:
                # Iterate over the interactions of the current type
                for interaction in interaction_elem:
                    resnr = interaction.find("resnr")
                    if resnr is not None:
                        resnr = resnr.text
                        restype = interaction.find("restype").text
                        data.append(
                            {
                                "resnr": resnr,
                                "restype": restype,
                                "interaction_type": interaction_type,
                            }
                        )
    # Create a DataFrame from the data list
    df = pd.DataFrame(data)
    # Remove duplicates to count interactions on a per residue basis
    df = df.drop_duplicates(keep="first")
    return df


def merge_dfs(df0, df1):
    """
    Merges the reference and query data frames and adds columns that state whether
    or not the residue is unique to the reference or query. This is what is used
    to calculate the Tanimoto coefficient.
    """
    # Rename the dataframes so that their info is distinct
    df0 = df0.rename(
        columns={
            "resnr": "resnr_1",
            "restype": "restype_1",
            "interaction_type": "interaction_type_1",
        }
    )
    df1 = df1.rename(
        columns={
            "resnr": "resnr_2",
            "restype": "restype_2",
            "interaction_type": "interaction_type_2",
        }
    )
    # Merge and add new columns that state T/F depending on if the row value is unique
    merged_df = pd.merge(
        df0,
        df1,
        left_on=["resnr_1", "restype_1", "interaction_type_1"],
        right_on=["resnr_2", "restype_2", "interaction_type_2"],
        how="outer",
    )
    merged_df = merged_df.assign(
        X=merged_df["resnr_1"].notnull(), Y=merged_df["resnr_2"].notnull()
    )
    merged_df["resnr_1"].update(merged_df["resnr_2"])
    merged_df["restype_1"].update(merged_df["restype_2"])
    merged_df["interaction_type_1"].update(merged_df["interaction_type_2"])
    # Drop extra columns and rename remaining columns
    merged_df = merged_df.drop(columns=["resnr_2", "restype_2", "interaction_type_2"])
    merged_df = merged_df.rename(
        columns={
            "resnr_1": "resnr",
            "restype_1": "restype",
            "interaction_type_1": "interaction",
            "X": "ref",
            "Y": "query",
        }
    )
    merged_df[["ref", "query"]] = merged_df[["ref", "query"]].astype(int)

    return merged_df


def similarity_plifs(plif_ref, plif_query):
    """
    Calcaultes the Tanimoto Similarity between the reference and query PLIFs.
    """
    sim = DataStructs.TanimotoSimilarity(plif_ref, plif_query)
    return sim


def get_vdw_contacts(pdb_file, ligand_name):
    """
    Calculates the van der Waals contacts between the protein and ligand.
    """
    # Create dictionary of vdw radii for common atoms
    vdw_radii = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}
    # Create empty lists for storing protein and ligand coordinates
    protein_coordinates = []
    ligand_coordinates = []
    # Create a temporary df for storing atom coordinates, distances, and vdw info
    df = pd.DataFrame(
        columns=[
            "HETATM",
            "LIG x",
            "LIG y",
            "LIG z",
            "ATOM",
            "restype",
            "resnr",
            "RES x",
            "RES y",
            "RES z",
            "DIST",
            "vdw_radii",
            "vdw_interaction",
        ]
    )
    # Read the atom lines of the pdb file to analyze
    with open(pdb_file) as f:
        lines = f.readlines()
        for line in lines:
            # Extract atom info
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_number = int(line[22:26].strip())
                res_name = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                # Save the protein and ligand atom info separately
                if res_name == ligand_name:
                    ligand_coordinates.append((atom_name, x, y, z))
                else:
                    protein_coordinates.append(
                        (atom_name, res_name, res_number, x, y, z)
                    )
    # Calculate euclidean distance using euclidean3d function
    # Set coordinate objects to enumerate so they can be iterated through
    for i, protein_row in enumerate(protein_coordinates):
        prot_x, prot_y, prot_z = protein_row[3:6]
        for j, ligand_row in enumerate(ligand_coordinates):
            lig_x, lig_y, lig_z = ligand_row[1:4]
            dist = euclidean3d([prot_x, prot_y, prot_z], [lig_x, lig_y, lig_z])
            # Threshold set to 4.5 angstroms to ensure maximum vdw_radii
            # are retained with additional atoms
            if dist < 4.5:
                temp_df = pd.DataFrame(
                    {
                        "HETATM": [ligand_row[0]],
                        "LIG x": [lig_x],
                        "LIG y": [lig_y],
                        "LIG z": [lig_z],
                        "ATOM": [protein_row[0]],
                        "restype": [protein_row[1]],
                        "resnr": [protein_row[2]],
                        "RES x": [prot_x],
                        "RES y": [prot_y],
                        "RES z": [prot_z],
                        "DIST": [dist],
                    }
                )
                df = pd.concat([df, temp_df], ignore_index=True)

    df["vdw_radii"] = df.apply(
        lambda row: vdw_radii.get(row["HETATM"][0], 0)
        + vdw_radii.get(row["ATOM"][0], 0)
        + 0.6,
        axis=1,
    )
    df["vdw_interaction"] = df.apply(
        lambda row: True if row["DIST"] < row["vdw_radii"] else False, axis=1
    )

    vdw_df = df[df["vdw_interaction"] == True]
    vdw_df = vdw_df.groupby(["resnr", "restype"]).first().reset_index()
    vdw_df = (
        vdw_df[["resnr", "restype"]]
        .rename(columns={"resnr": "resnr", "restype": "restype"})
        .assign(interaction_type="vdw_contact")
    )
    return vdw_df


def sim_matrix(tanimoto_dict):
    """
    Generates a similarity matrix from a dictionary of Tanimoto similarities.
    """
    # Split the dictionary keys into sets
    keys = list(
        set(
            [i.split(" & ")[0] for i in tanimoto_dict.keys()]
            + [i.split(" & ")[1] for i in tanimoto_dict.keys()]
        )
    )
    # Create an empty matrix
    matrix = np.zeros((len(keys), len(keys)))
    # Fill the matrix with Tanimoto similarities
    for k in tanimoto_dict.keys():
        pdb1 = k.split(" & ")[0]
        pdb2 = k.split(" & ")[1]
        matrix[list(keys).index(pdb1), list(keys).index(pdb2)] = tanimoto_dict[k]
    # Create a dataframe from the matrix
    df = pd.DataFrame(matrix, index=keys, columns=keys)
    # Order the dataframe alphabetically and starting with numbers
    df = df.reindex(sorted(df.columns, key=lambda x: (x.isdigit(), x)), axis=1)
    df = df.reindex(sorted(df.index, key=lambda x: (x.isdigit(), x)), axis=0)
    # Write the data frame to a CSV file named "tanimoto_similarity_matrix.csv"
    df.to_csv("tanimoto_similarity_matrix.csv")
    # Create a heatmap
    sns.heatmap(df, cmap="inferno", cbar_kws={"label": "Tanimoto Similarity"})
    # Adjust the x axis labels
    plt.xticks(rotation=45, va="center", ha="right", rotation_mode="anchor")
    # Reduce the size of the axis labels
    plt.tick_params(axis="both", which="major", labelsize=8)
    # Make sure the x and y axis labels are visible
    plt.tight_layout()
    # Save the plot as a png file
    plt.savefig("tanimoto_similarity_matrix.png", dpi=300)
    # Display the heatmap
    plt.show()


def check_exit(input_str):
    """
    Checks if the user inputted "exit" and return True if they did.
    This allows the user to exit the program at any prompt.
    """
    if isinstance(input_str, str):
        input_str = input_str.lower()
        if input_str == "exit":
            return True
        return False


def PLIP_PLIFs():
    """
    Main function
    """
    print(
        '\nNOTE: All of the PLIP XML reports and protonated PDB files must be in the same directory.\nThe PDB files must all end with "_protonated.pdb".\nType "exit" at any prompt to exit the program.\n'
    )
    # Ask the user to enter the working directory
    pwd = input('Enter your working directory (type "cd" or "current" to skip): ')
    # IF the user enters 'current' or 'cd', do not change the working directory
    if pwd == "current" or pwd == "cd":
        pwd = os.getcwd()
    if check_exit(pwd):
        return
    # Keep asking the user to enter the working directory until it exists
    while not os.path.exists(pwd):
        print(
            "Error: This directory does not exist. Make sure that the directory path is correct.\n"
        )
        pwd = input("Enter your working directory: ")
        if check_exit(pwd):
            return
    # Change the working directory
    os.chdir(pwd)
    # Ask the user to enter the ligand name
    ligand = input("Enter the ligand ID as it is found within the PDB files: ")
    # Ask the user to enter the name of the reference PDB file
    ref_pdb = input("Enter the name of the protonated reference PDB file: ")
    if check_exit(ref_pdb):
        return
    if not ref_pdb.endswith(".pdb"):
        ref_pdb += ".pdb"
    # Check if the ref_pdb ends with "_protonated.pdb"
    while not ref_pdb.endswith("_protonated.pdb"):
        print(
            'Error: The reference PDB file name must be protonated and the file name must end with "_protonated.pdb".'
        )
        ref_pdb = input("Enter the file name of the protonated reference PDB file: ")
        if not ref_pdb.endswith(".pdb"):
            ref_pdb += ".pdb"
        if check_exit(ref_pdb):
            return
    # Check if the file exists
    while not os.path.exists(ref_pdb):
        print(
            'Error: This file does not appear to exist. Make sure that the file name ends with "_protonated.pdb". The entry is case sensitive.'
        )
        ref_pdb = input("Enter the file name of the protonated reference PDB file: ")
        if not ref_pdb.endswith(".pdb"):
            ref_pdb += ".pdb"
        if check_exit(ref_pdb):
            return
    # Get the base name of the reference PDB file without the '_protonated' ending
    ref_pdb_name = ref_pdb.split("_protonated")[0].split(".")[0]
    # Create dictionaries to store the PLIFS, merged data frames, and Tanimoto coefficients
    plifs = {}
    tanimoto_dict = {}
    interaction_dfs = {}
    # Start a timer
    start_time = time.time()
    # Loop through all the PDB and XML files in the working directory
    for file in os.listdir(pwd):
        if file.endswith("_protonated.pdb"):
            pdb = file
            # Get the base name of the PDB file without the '_protonated' ending
            pdb_name = pdb.split("_protonated")[0].split(".")[0]
            # Define the PLIP XML report file name as the PDB file name with the extension changed to .xml and without the '_protonated' ending
            xml = pdb_name + ".xml"
            print("Parsing PLIP XML report for " + pdb_name + "...")
            # Parse the PLIP XML report
            plip_df = parse_plip_xml(xml)
            # Get the Van der Waals contacts for the ligand
            print("Calculating Van der Waals contacts for " + pdb_name + "...")
            vdw_df = get_vdw_contacts(pdb, ligand)
            # Combine the PLIP and Van der Waals data frames
            df0 = pd.concat([plip_df, vdw_df])
            df0 = df0.reset_index(drop=True)
            # Add the df to the interactions dictionary
            interaction_dfs[pdb_name] = df0
    # Loop through the plifs dictionary and with all possible combinations of PLIFs, merge the data frames
    for k1, v1 in interaction_dfs.items():
        for k2, v2 in interaction_dfs.items():
            # Merge the PLIFs
            df1 = merge_dfs(v1, v2)
            # Add the merged data frame to the plifs dictionary
            plifs[k1 + " & " + k2] = df1
    # Write only the merged data frames containing the reference PLIF to a text file
    for k, v in plifs.items():
        if ref_pdb_name in k:
            with open(k + "_PLIF.txt", "w") as f:
                f.write(v.to_string())
    # Delete the text files that end with ref_pdb_name + '_PLIF.txt'
    for file in os.listdir(pwd):
        if file.endswith(ref_pdb_name + "_PLIF.txt"):
            os.remove(file)
    # Change the titles of the text files so that everything before " & " is removed
    for file in os.listdir(pwd):
        if file.endswith("_PLIF.txt"):
            new_file = file.split(" & ")[1]
            os.rename(file, new_file)
    # Loop through all of the merged data frames and calculate the Tanimoto coefficient
    print("Calculating Tanimoto coefficients...")
    for k, v in plifs.items():
        # Convert the columns to bit strings
        bit_string1 = str(v["ref"].to_numpy())
        bit_string2 = str(v["query"].to_numpy())
        plif1_bit = rdkit.DataStructs.cDataStructs.CreateFromBitString(bit_string1)
        plif2_bit = rdkit.DataStructs.cDataStructs.CreateFromBitString(bit_string2)
        # Calculate the Tanimoto coefficient
        tanimoto = similarity_plifs(plif1_bit, plif2_bit)
        tanimoto = round(tanimoto, 3)
        # Add the Tanimoto coefficient to the tanimoto dictionary
        if k + " & " + k not in tanimoto_dict:
            tanimoto_dict[k + " & " + k] = tanimoto
        else:
            pass
    print("\nWorking on generating a similarity matrix...")
    # Run the function to generate a similarity matrix
    sim_matrix(tanimoto_dict)
    print("Moving output files to PLIF_files folder...")
    # Create a new folder in the working directory to store the non-structure files
    if not os.path.exists(os.path.join(pwd, "PLIF_files")):
        os.makedirs(os.path.join(pwd, "PLIF_files"))
    # Move the text files and XML files to the new folder
    try:
        for file in os.listdir(pwd):
            if file.endswith("PLIF.txt") or "tanimoto_similarity_matrix" in file:
                shutil.move(file, os.path.join(pwd, "PLIF_files"))
        # Inform the user where the files were moved to
        print(
            "All PLIF files moved to the PLIF_files folder in the current working directory."
        )
    except:
        print(
            "Identical files were found in the output folder destination. No files were moved."
        )
    # End the timer
    end_time = time.time()
    # Calculate the total time
    total_time = end_time - start_time
    # Print the total time in minutes
    print("\nTotal run time: " + str(round(total_time / 60, 2)) + " minutes")


if __name__ == "__main__":
    PLIP_PLIFs()
