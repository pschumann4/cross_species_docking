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
    pdb_dir = input("Enter the directory containing the PDB models: ")

    # Check that the working directory exists
    while not os.path.exists(pdb_dir):
        print("The working directory does not exist.")
        pdb_dir = input('Enter the directory (or type "exit"): ')

    # Change the working directory
    os.chdir(pdb_dir)

    # Search up to 3 parent levels for a 'details' subdirectory
    output_dir = pdb_dir  # default fallback
    search_path = pdb_dir
    for level in range(4):  # 0 = current dir, 1-3 = parent levels
        details_candidate = os.path.join(search_path, "details")
        if os.path.isdir(details_candidate):
            output_dir = details_candidate
            print(f"Found 'details' directory at: {details_candidate}")
            break
        search_path = os.path.dirname(search_path)
    else:
        print("No 'details' directory found within 3 parent levels. Saving to working directory.")

    # Search the directory for a file that starts with "ref_"
    ref_pdb = ""
    pdb_files = [i for i in os.listdir(pdb_dir) if i.endswith(".pdb")]
    ref = "n"
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

    files = [i for i in os.listdir(pdb_dir) if i.endswith(".pdb") and "model" in i and i != ref_pdb]

    # Calculate RMSDs between reference ligand and all query ligands
    rmsd_dict = {}
    files = [i for i in os.listdir(pdb_dir) if i.endswith(".pdb") and "model" in i and i != ref_pdb]
    
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
    with open(os.path.join(output_dir, "ligand_rmsd.txt"), "w") as f:
        # Add a title line
        f.write("Ligand RMSD values\n" + "\n")
        # Add the RMSD values for each PDB file
        for pdb2, rmsd in rmsd_dict.items():
            f.write("%s: %s\n" % (pdb2, rmsd))

    # Remove the quotes from the text file and the spaces between the residue name and commas
    rmsd_txt_path = os.path.join(output_dir, "ligand_rmsd.txt")
    with open(rmsd_txt_path, "r") as f:
        lines = f.readlines()
    with open(os.path.join(output_dir, "ligand_rmsd.txt"), "w") as f:
        for line in lines:
            line = line.replace("'", "")
            line = line.replace(", ", ",")
            f.write(line)

    # Check that the text file was created
    if os.path.exists(rmsd_txt_path):
        print(
            "The ligand RMSD values were written to {}.".format(rmsd_txt_path)
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
    plt.savefig(os.path.join(output_dir, "ligand_rmsd.png"), dpi=300)
    plt.show()

if __name__ == "__main__":
    ligand_rmsd()
