"""
 This script will calculate the average pose of a ligand from a set of PDB files.  
 The average pose will be written to a PDB file.                                  
"""


import os

import numpy as np


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


def generate_ref_pose():
    """
    This is the main function. It calculates the average pose of a ligand 
    from a set of PDB files.
    """
    # Ask the user for the working directory
    pwd = input("Enter the directory containing the PDB models: ")
    if check_exit(pwd):
        return
    while not os.path.exists(pwd):
        print("The working directory does not exist.")
        pwd = input('Enter the directory (or type "exit"): ')
        if check_exit(pwd):
            return
    os.chdir(pwd)
    # Ask the user for the ligand name
    ligand = input("Enter the ligand ID as it is found within the PDB files: ")
    if check_exit(ligand):
        return
    # Loop through all of the PDB files in the working directory
    # and check that they contain the ligand
    pdb_files = []
    no_lig = []
    for file in os.listdir():
        if file.endswith(".pdb") and "model" in file:
            with open(file, "r") as f:
                for line in f:
                    if line.startswith("HETATM") and line[17:20].strip() == ligand:
                        pdb_files.append(file)
                        break
                    if line.startswith("HETATM") and line[17:20].strip() != ligand:
                        if file not in no_lig:
                            no_lig.append(os.path.basename(file))
    # Print the names of the PDB files that do not contain the ligand
    if len(no_lig) > 0:
        print("WARNING: The following PDB files do not contain the ligand: ")
        for file in no_lig:
            print(file)
    ligand_coords = []
    # Iterate through the PDB files to get ligand coordinates
    for file in pdb_files:
        coords = []
        atom_type = []
        with open(file, "r") as f:
            for line in f:
                if line.startswith("HETATM") and line[17:20].strip() == ligand:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom_type.append(line[12:16].strip())
                    coords.append([x, y, z])
            ligand_coords.append(coords)
    # Convert the coordinates to a numpy arrays
    ligand_coords = [np.array(coords) for coords in ligand_coords]
    # Choose model to serve as a reference (to maintain the original atom order)
    ref = ligand_coords[0]
    # Calculate the average pose by taking the mean of the coordinates
    avg_pose = np.mean(ligand_coords, axis=0)
    # Update the coordinates of the ref molecule with the average pose
    ref = avg_pose + (ref - avg_pose)
    # Write the coordinates to a pdb file
    with open("avg_pose.pdb", "w") as f:
        for i, coords in enumerate(ref):
            f.write(
                f"HETATM{i+1:>5} {atom_type[i]:<4} {ligand:>3} A{1:>4}    {coords[0]:>8.3f}{coords[1]:>8.3f}{coords[2]:>8.3f}\n"
            )
        f.write("END")
    print("The average pose has been written to avg_pose.pdb")

if __name__ == "__main__":
    generate_ref_pose()
