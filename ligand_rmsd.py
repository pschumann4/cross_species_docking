"""
 This script will calculate the RMSD between the ligand in a reference PDB file and the ligand in 
 a set of other PDB files.                                                                        
 The output will be a text file containing the RMSD values for each PDB file.                     
"""
import os
import shutil

import numpy as np
import pandas as pd


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


def ligand_rmsd():
    """
    Main function
    """
    # Prompt user for the working directory
    pwd = input("Enter the directory containing the PDB models: ")
    if check_exit(pwd):
        return

    # Check that the working directory exists
    while not os.path.exists(pwd):
        print("The working directory does not exist.")
        pwd = input('Enter the directory (or type "exit"): ')
        if check_exit(pwd):
            return

    # Change the working directory
    os.chdir(pwd)

    # Prompt user for the reference PDB file
    ref = input("Enter the name of the reference PDB file: ")
    if check_exit(ref):
        return

    # Check if the PDB file ends with '.pdb'
    if not ref.endswith(".pdb"):
        ref += ".pdb"

    # Check that the PDB file exists
    while not os.path.exists(ref):
        print(
            "The pdb file you entered does not appear to exist. The entry is case sensitive."
        )
        ref = input('Enter the name of the reference PDB file (or type "exit"): ')
        if not ref.endswith(".pdb"):
            ref += ".pdb"
        if check_exit(ref):
            return

    # Prompt user for the ligand name
    ligand = input("Enter the ligand ID as it is found within the PDB files: ")
    if check_exit(ligand):
        return

    # Read the PDB file and store the HETATM lines for the specified ligand in a list
    hetatm1 = []
    with open(ref, "r") as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == ligand:
                hetatm1.append(line)
    # Check that the ligand is present in the PDB file
    if len(hetatm1) == 0:
        print("Error: The ligand was not found in the PDB file!")
        return
    # Create a list to store the coordinates of the ligand
    coords1 = []
    atoms1 = []

    # Iterate through the HETATM lines
    for line in hetatm1:
        # Extract the x, y, and z coordinates from the line (unless it's a hydrogen)
        if "H" not in line[12:16]:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            # Extract atom name
            atoms1.append(line[12:16].strip())
            # Add coordinates to the list
            coords1.append([x, y, z])
    # Convert the list to a NumPy array
    coords1 = np.array(coords1)
    # Create a dictionary to store the RMSD values using the base name of pdb2 as the key
    rmsd_dict = {}

    # Calculate the RMSDs between the reference ligand and all query ligands
    for file in os.listdir(pwd):
        # Check if the file is a PDB file containing the word "model" and is not the reference PDB file
        if file.endswith(".pdb") and "model" in file and file != ref:
            # Get the base name of the PDB file
            base = os.path.basename(file)
            base = os.path.splitext(base)[0]
            # Read the PDB file and store the HETATM lines for the specified ligand in a list
            hetatm2 = []
            with open(file, "r") as f:
                for line in f:
                    if line.startswith("HETATM") and line[17:20].strip() == ligand:
                        hetatm2.append(line)
            # Check that the ligand is present in the PDB file, if not, skip the file
            if len(hetatm2) == 0:
                print("The ligand was not found in " + file + ".")
            else:
                # Create a list to store the coordinates of the ligand
                coords2 = []
                atoms2 = []
                # Iterate through the HETATM lines
                for line in hetatm2:
                    # Extract the x, y, and z coordinates from the line (unless it's a hydrogen)
                    if "H" not in line[12:16]:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        # Extract atom name
                        atoms2.append(line[12:16].strip())
                        # Add the coordinates to the list
                        coords2.append([x, y, z])
                # Convert the list to a NumPy array
                coords2 = np.array(coords2)
                # Check that the atoms in the two PDB files are the same
                if atoms1 != atoms2:
                    print(
                        "ERROR: The atoms in "
                        + ref
                        + " and "
                        + file
                        + " are not identical.\nThis must be corrected before the RMSD can be calculated properly."
                    )
                    # Print the lists of atoms in each PDB file
                    print("Atoms in " + ref + ": " + str(atoms1))
                    print("Atoms in " + file + ": " + str(atoms2))
                    return
                # Perform RMSD calculation
                diff = coords1 - coords2
                rmsd = np.sqrt(np.sum(diff**2) / len(coords1))
                # Round the rmsd value to 3 decimal places
                rmsd = round(rmsd, 3)
                # Add the RMSD value to the dictionary
                rmsd_dict[base] = rmsd

    # Write the RMSD values to a file
    with open("ligand_rmsd.txt", "w") as f:
        # Add a title line
        f.write("Ligand RMSD values\n" + "\n")
        # Add the RMSD values for each PDB file
        for pdb2, rmsd in rmsd_dict.items():
            f.write("%s: %s\n" % (pdb2, rmsd))

    # Remove the quotes from the text file and the spaces between the residue name and commas
    with open("ligand_rmsd.txt", "r") as f:
        lines = f.readlines()
    with open("ligand_rmsd.txt", "w") as f:
        for line in lines:
            line = line.replace("'", "")
            line = line.replace(", ", ",")
            f.write(line)

    # Check that the text file was created
    if os.path.exists("ligand_rmsd.txt"):
        print(
            "\nThe ligand RMSD values were written to ligand_rmsd.txt in your working directory."
        )

    # Determine best and worst poses?
    filter_poses = input(
        "\nWould you like to sort the models as " "best/worst based on RMSD? (y/n): "
    )
    filter_poses = filter_poses.lower()
    while filter_poses not in ["y", "n"]:
        filter_poses = input("Please enter y or n: ")
        filter_poses = filter_poses.lower()
    if filter_poses == "y":
        lig_rmsd_file = "ligand_rmsd.txt"
        # Initialize a df to hold the ligand RMSD values, species, and model
        lig_rmsd_df = pd.DataFrame(columns=["lig_rmsd", "species", "model"])
        # Read the ligand RMSD file
        with open(lig_rmsd_file, "r") as f:
            lines = f.readlines()
            lines = [line.split() for line in lines]
            lig_rmsds = [line[1] for line in lines[2:]]
            species = [line[0].split("_")[0] for line in lines[2:]]
            # Find the model number by selecting the number after the word "model"
            model = [line[0].split("model")[1].split(":")[0] for line in lines[2:]]
        # Add the ligand RMSD values to the dataframe
        lig_rmsd_df["lig_rmsd"] = lig_rmsds
        lig_rmsd_df["species"] = species
        lig_rmsd_df["model"] = model
        # For each species in the dataframe, find the best and worst pose
        best_poses = []
        worst_poses = []
        for species in lig_rmsd_df["species"].unique():
            species_df = lig_rmsd_df[lig_rmsd_df["species"] == species]
            best_pose = species_df[
                species_df["lig_rmsd"] == species_df["lig_rmsd"].min()
            ]
            worst_pose = species_df[
                species_df["lig_rmsd"] == species_df["lig_rmsd"].max()
            ]
            best_poses.append(best_pose)
            worst_poses.append(worst_pose)
        # Concatenate the best and worst poses into a single dataframe
        best_poses = pd.concat(best_poses)
        worst_poses = pd.concat(worst_poses)
        # Remove the index from the best and worst poses df
        best_poses.reset_index(drop=True, inplace=True)
        worst_poses.reset_index(drop=True, inplace=True)
        # Convert the lig_rmsd column to a float
        best_poses["lig_rmsd"] = best_poses["lig_rmsd"].astype(float)
        worst_poses["lig_rmsd"] = worst_poses["lig_rmsd"].astype(float)
        # Convert the pose column to an integer
        best_poses["model"] = best_poses["model"].astype(int)
        worst_poses["model"] = worst_poses["model"].astype(int)

        # Write the best and worst poses to a XLSX file
        with pd.ExcelWriter("best_and_worst_poses.xlsx") as writer:
            best_poses.to_excel(writer, sheet_name="best_poses")
            worst_poses.to_excel(writer, sheet_name="worst_poses")
        print(
            "\nThe best and worst poses have been written "
            'to "best_and_worst_poses.xlsx".'
        )

        # Create a folder called "filtered_models" to hold the best and worst poses
        filtered_models_dir = os.path.join(pwd, "filtered_models")
        if not os.path.exists(filtered_models_dir):
            os.mkdir(filtered_models_dir)

        # Copy the best and worst poses to the "filtered_models" folder
        for model in os.listdir(pwd):
            if model.endswith(".pdb") and "model" in model:
                # Get the species and model number from the model name
                species = model.split("_")[0]
                model_num = model.split("model")[1].split(".pdb")[0]
                for pose in best_poses.itertuples():
                    if pose.species == species and pose.model == int(model_num):
                        shutil.copy(model, filtered_models_dir)
                for pose in worst_poses.itertuples():
                    if pose.species == species and pose.model == int(model_num):
                        shutil.copy(model, filtered_models_dir)

        print(
            '\nDONE! The best and worst poses have been copied to the "filtered_models" folder.'
        )


if __name__ == "__main__":
    ligand_rmsd()
