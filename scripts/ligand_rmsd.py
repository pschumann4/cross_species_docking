import os
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def ligand_rmsd():
    """
    Calculate the RMSD between a reference ligand and ligands in other PDB files.
    For this function to work properly, the bound ligand in each model must be
    named the same as the reference ligand.

    Input:
    - PDB files containing the ligands to be compared.
    - The reference PDB file containing the reference ligand.
    - The ligand ID as it is found within the PDB files.

    Output:
    - A text file containing the RMSD values for each PDB file.
    - A histogram showing the distribution of RMSD values.
    - Optionally, the best and worst poses can be filtered and copied to a new directory.
    """
    # Prompt user for the working directory
    cwd = input("Enter the directory containing the PDB models: ")

    # Check that the working directory exists
    while not os.path.exists(cwd):
        print("The working directory does not exist.")
        cwd = input('Enter the directory (or type "exit"): ')

    # Change the working directory
    os.chdir(cwd)

    # Search the directory for a file that starts with "ref_"
    ref_pdb = ""
    pdb_files = [i for i in os.listdir(cwd) if i.endswith(".pdb")]
    for file in pdb_files:
        if file.startswith("ref_"):
            ref = input("Is {} the reference PDB file? (y/n): ".format(file))
            ref = ref.lower()
            # check for valid input
            while ref.lower() not in ["y", "n"]:
                ref = input("Please enter y or n: ")
            if ref == "y":
                ref_pdb = file
                break
    if ref == "n":
        ref_pdb = input("Enter the name of the reference PDB file (if none, type 'random'): ")

        if not ref_pdb.endswith(".pdb"):
            ref_pdb += ".pdb"
    
    # Check that the reference file exists and if not ask for it again
    while not os.path.exists(ref_pdb):
        print("The reference file does not exist.")
        ref_pdb = input("Enter the name of the reference PDB file: ")

        if not ref_pdb.endswith(".pdb"):
            ref_pdb += ".pdb"

    # Prompt user for the ligand name
    ligand = input("Enter the ligand ID as it is found within the PDB files: ")

    # Read the PDB file and store the HETATM lines for the specified ligand in a list
    hetatm1 = []
    with open(ref_pdb, "r") as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == ligand:
                hetatm1.append(line)
    # Check that the ligand is present in the PDB file
    while not hetatm1:
        print("ERROR: The ligand entered was not found in the PDB file.")
        # List all the ligands found in the PDB file
        ligands = set()
        with open(ref_pdb, "r") as f:
            for line in f:
                if line.startswith("HETATM"):
                    ligands.add(line[17:20].strip())
        print("The following ligands were found in the PDB file:")
        print(ligands)
        ligand = input("Please select one of the ligands from the list: ")
        hetatm1 = []
        with open(ref_pdb, "r") as f:
            for line in f:
                if line.startswith("HETATM") and line[17:20].strip() == ligand:
                    hetatm1.append(line)

    # Create a list to store the coordinates of the ligand
    coords1 = []
    elements1 = []

    # Iterate through the HETATM lines
    for line in hetatm1:
        # Extract the x, y, and z coordinates from the line
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        # Extract atom name
        atom_name = line[12:16].strip()
        # Extract element from atom name
        element = atom_name[:1]
        if element == "H":
            continue
        elements1.append(element)
        # Add coordinates to the list
        coords1.append([x, y, z])

    # Convert the list to a NumPy array
    coords1 = np.array(coords1)

    # Sort the coordinates
    coords1 = np.sort(coords1, axis=0)

    # Calculate RMSDs between reference ligand and all query ligands
    rmsd_dict = {}

    files = [i for i in os.listdir(cwd) if i.endswith(".pdb") and "model" in i and i != ref_pdb]

    # Calculate RMSDs between reference ligand and all query ligands
    rmsd_dict = {}
    files = [i for i in os.listdir(cwd) if i.endswith(".pdb") and "model" in i and i != ref_pdb]
    
    for file in files:
        base = os.path.basename(file)
        base = os.path.splitext(base)[0]
        hetatm2 = []
        
        with open(file, "r") as f:
            for line in f:
                if line.startswith("HETATM") and line[17:20].strip() == ligand:
                    hetatm2.append(line)
                    
        if len(hetatm2) == 0:
            print("The ligand was not found in " + file + ".")
            continue
            
        coords2 = []
        elements2 = []
        
        for line in hetatm2:
            if "H" not in line[12:16]:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_name = line[12:16].strip()
                element = atom_name[:1]
                if element == "H":
                    continue
                elements2.append(element)
                coords2.append([x, y, z])
                
        coords2 = np.array(coords2)
        coords2 = np.sort(coords2, axis=0)
        
        # Check that the number of atoms and element types match
        if len(elements1) != len(elements2):
            print(f"ERROR: Different number of atoms in {ref_pdb} ({len(elements1)}) and {file} ({len(elements2)})")
            continue
            
        # Sort elements and check if they match
        if sorted(elements1) != sorted(elements2):
            print(f"ERROR: Different element types in {ref_pdb} and {file}")
            print(f"Elements in {ref_pdb}: {sorted(elements1)}")
            print(f"Elements in {file}: {sorted(elements2)}")
            continue
            
        # Perform RMSD calculation
        diff = coords1 - coords2
        rmsd = np.sqrt(np.sum(diff**2) / len(coords1))
        rmsd_dict[base] = round(rmsd, 3)

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
            "The ligand RMSD values were written to ligand_rmsd.txt "
            "in your working directory."
        )

    # Plot the RMSD values as a histogram
    rmsd_values = []
    # Iterate through the dictionary and add the RMSD values to the list
    for rmsd in rmsd_dict.values():
        rmsd_values.append(rmsd)
    plt.hist(rmsd_values, bins=30)
    plt.xlabel("RMSD (Å)", fontsize=18)
    plt.ylabel("Frequency", fontsize=18)
    for patch in plt.gca().patches:
        patch.set_facecolor("#A4D4F7")
        patch.set_edgecolor("black")
    plt.tight_layout()
    plt.savefig("ligand_rmsd.png", dpi=300)
    plt.show()

    # Determine best and worst poses?
    filter_poses = input(
        "\nWould you like to filter the best/worst models based on RMSD? (y/n): "
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
            model = [line[0].split("model")[1].split(":")[0] for line in lines[2:]]
        # Add the ligand RMSD values to the dataframe
        lig_rmsd_df["lig_rmsd"] = lig_rmsds
        lig_rmsd_df["lig_rmsd"] = lig_rmsd_df["lig_rmsd"].astype(float)
        lig_rmsd_df["species"] = species
        lig_rmsd_df["model"] = model

        # Rmove any models that have a ligand RMSD > 10 Å
        lig_rmsd_df = lig_rmsd_df[lig_rmsd_df["lig_rmsd"] < 10]

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
            "The best and worst poses have been written "
            'to "best_and_worst_poses.xlsx".'
        )

        # Create a folder called "filtered_models" to hold the best and worst poses
        filtered_models_dir = os.path.join(cwd, "filtered_models")
        if not os.path.exists(filtered_models_dir):
            os.mkdir(filtered_models_dir)

        # Copy the best and worst poses to the "filtered_models" folder
        for model in os.listdir(cwd):
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
            if model.startswith("ref_"):
                shutil.copy(model, filtered_models_dir)

        print(
            "DONE! The best and worst poses have been copied "
            'to the "filtered_models" folder.'
        )

        # Edit the ligand_rmsd.txt file so that it only includes the best and worst poses
        with open("ligand_rmsd.txt", "r") as f:
            lines = f.readlines()
        written_lines = set()
        with open("filtered_ligand_rmsd.txt", "w") as f:
            # Write the first two lines of the file
            for line in lines[:2]:
                f.write(line)
            
            # Write the best and worst poses to the file
            for line in lines[2:]:
                # Create a tuple of the identifying parts of the line to use as a unique key
                species = line.split("_")[0]
                model = line.split("model")[1].split(":")[0]
                line_key = (species, model)
                
                # Check if this line should be included (in best or worst poses)
                should_include = False
                
                for pose in best_poses.itertuples():
                    if pose.species == species and str(pose.model) == model:
                        should_include = True
                        break
                        
                if not should_include:
                    for pose in worst_poses.itertuples():
                        if pose.species == species and str(pose.model) == model:
                            should_include = True
                            break
                
                # If this line should be included and hasn't been written yet, write it
                if should_include and line_key not in written_lines:
                    f.write(line)
                    written_lines.add(line_key)

        print(
            "A filtered ligand RMSD file has been created to only include "
            "the best and worst poses."
        )

        # Re-create the histogram with the best and worst poses highlighted
        rmsd_values = []
        # Iterate through the dictionary and add the RMSD values to the list
        for rmsd in rmsd_dict.values():
            rmsd_values.append(rmsd)
        plt.hist(rmsd_values, bins=30)
        plt.xlabel("RMSD (Å)", fontsize=18)
        plt.ylabel("Frequency", fontsize=18)
        # Add vertical lines to show the range of the best and worst poses
        max_best = best_poses["lig_rmsd"].max()
        min_best = best_poses["lig_rmsd"].min()
        plt.axvline(
            x=max_best,
            color="#EB0744",
            linestyle="--",
            label="Best poses",
            linewidth=2,
        )
        plt.axvline(
            x=min_best,
            color="#EB0744",
            linestyle="--",
            linewidth=2,
        )
        max_worst = worst_poses["lig_rmsd"].max()
        min_worst = worst_poses["lig_rmsd"].min()
        plt.axvline(
            x=max_worst,
            color="#062576",
            linestyle="--",
            label="Worst poses",
            linewidth=2,
        )
        plt.axvline(
            x=min_worst,
            color="#062576",
            linestyle="--",
            linewidth=2,
        )
        for patch in plt.gca().patches:
            patch.set_facecolor("#A4D4F7")
            patch.set_edgecolor("black")
        plt.legend()
        plt.tight_layout()
        plt.savefig("labeled_ligand_rmsd.png", dpi=300)
        plt.show()

        filter_vina = input(
            "\nWould you like to filter the Vina output according to RMSD "
            "as well? (y/n): "
        )
        filter_vina = filter_vina.lower()
        while filter_vina not in ["y", "n"]:
            filter_vina = input("Please enter y or n: ")
            filter_vina = filter_vina.lower()
        if filter_vina == "y":
            vina_logs = input("Enter the path to the Vina output logs: ")
            while not os.path.exists(vina_logs):
                vina_logs = input("Please enter a valid path: ")
            os.chdir(vina_logs)
            filtered_vina_dir = os.path.join(vina_logs, "filtered_vina_output")
            if not os.path.exists(filtered_vina_dir):
                os.mkdir(filtered_vina_dir)
            vina_files = [i for i in os.listdir(vina_logs) if i.endswith(".pdbqt") and "bound" in i]
            for file in vina_files:
                species = file.split("_")[0]
                if "ligand" in file:
                    model_num = file.split("ligand_")[1].split(".pdbqt")[0]
                if "flex" in file:
                    model_num = file.split("flex_")[1].split(".pdbqt")[0]
                for pose in best_poses.itertuples():
                    if pose.species == species and pose.model == int(model_num):
                        shutil.copy(file, filtered_vina_dir)
                for pose in worst_poses.itertuples():
                    if pose.species == species and pose.model == int(model_num):
                        shutil.copy(file, filtered_vina_dir)

        print(
            "DONE! The best and worst poses have been copied to "
            'the "filtered_vina_output" folder.'
        )


if __name__ == "__main__":
    ligand_rmsd()
