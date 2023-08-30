"""
 This function will essentially just replace the residue positions of the target  
 with the reference given that the specified residues are identical. However, it  
 will also adjust the atom positions by calculating the difference between the    
 reference and target backbone atoms and then adding that difference to the       
 target atoms.    
"""

import os
import string

import numpy as np


def modify_res_pos():
    """
    Main function
    """
    # Prompt user for the working directory
    wd = input("Enter the working directory: ")
    # Check that the directory exists
    while not os.path.exists(wd):
        print("The working directory does not exist.")
        wd = input("Enter the working directory: ")
    # Update the working directory
    os.chdir(wd)

    # Prompt user for the reference PDB file
    pdb_ref = input("Enter the name of the reference PDB file: ")
    if not pdb_ref.endswith(".pdb"):
        pdb_ref += ".pdb"

    # Read in the reference PDB file
    with open(pdb_ref, "r") as f:
        pdb_ref = f.readlines()

    # Prompt user for the target PDB file
    while True:
        pdb_tar = input("Enter the name of the target PDB file: ")
        if not pdb_tar.endswith(".pdb"):
            pdb_tar += ".pdb"
        # Get the name of the target PDB file by removing the extension
        pdb_tar_name = pdb_tar.split(".")[0]

        # Prompt user for the residues to be modified
        residues = input(
            "Enter a comma separated list of the residue numbers to be "
            "modified (with or without residue IDs): "
        )
        # Convert the comma-separated list of res to a list of res as ints
        if ", " in residues:
            residues = residues.split(", ")
        elif "," in residues:
            residues = residues.split(",")
        for i in range(len(residues)):
            if not residues[i].isdigit():
                residues[i] = residues[i].strip(string.ascii_letters)
        residues = [int(residue) for residue in residues]
        # Create a list of the new lines for the PDB file
        new_lines = []
        new_pdb = []

        # Read in the PDB files
        with open(pdb_tar, "r") as f:
            pdb_tar = f.readlines()
        for res in residues:
            ref_coords = {}
            # Extract the coordinates of the atoms in the ref PDB
            for line in pdb_ref:
                if line.startswith("ATOM"):
                    if int(line[22:26]) == res:
                        atom = line[12:16].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ref_coords[atom] = np.array([x, y, z])

            # Extract the coordinates of the atoms in the target PDB
            tar_coords = {}
            for line in pdb_tar:
                if line.startswith("ATOM"):
                    if int(line[22:26]) == res:
                        atom = line[12:16].strip()
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        tar_coords[atom] = np.array([x, y, z])

            # Calculate the centroids of the backbone atoms in the 
            # target and reference PDBs
            ref_centroid = np.mean(
                [ref_coords["N"], ref_coords["CA"], ref_coords["C"], ref_coords["O"]],
                axis=0,
            )
            tar_centroid = np.mean(
                [tar_coords["N"], tar_coords["CA"], tar_coords["C"], tar_coords["O"]],
                axis=0,
            )
            # Calculate the difference between the centroids
            diff = tar_centroid - ref_centroid

            # Adjust the residue coordinates by the difference between the centroids
            for atom in ref_coords:
                ref_coords[atom] += diff.reshape(ref_coords[atom].shape)
                ref_coords[atom] = ref_coords[atom].tolist()

            # Round the coordinates to 3 decimal places
            for atom in ref_coords:
                for i in range(len(ref_coords[atom])):
                    ref_coords[atom][i] = round(ref_coords[atom][i], 3)
                    ref_coords[atom][i] = str(ref_coords[atom][i])
                    # Check that each value has 3 decimal places, if not, 
                    # add a 0 to the end
                    if len(ref_coords[atom][i].split(".")[1]) < 3:
                        ref_coords[atom][i] += "0"
                    # Modify the strings so that they are 8 characters long
                    while len(ref_coords[atom][i]) < 8:
                        ref_coords[atom][i] = " " + ref_coords[atom][i]

                # Replace the coordinates of the atoms in the target PDB with 
                # the coordinates from the ref_coords dictionary
                for line in pdb_tar:
                    if line.startswith("ATOM"):
                        if int(line[22:26]) == res and line[12:16].strip() == atom:
                            # Replace the coordinates of that line with the 
                            # coordinates from the ref_coords dictionary
                            line = (
                                line[:30]
                                + "{:>8}".format(ref_coords[atom][0])
                                + "{:>8}".format(ref_coords[atom][1])
                                + "{:>8}".format(ref_coords[atom][2])
                                + line[54:]
                            )
                            new_lines.append(line)

        # Overwrite the corresponding lines in the target PDB with the 
        # lines in the new_pdb list
        for line in pdb_tar:
            if line.startswith("ATOM"):
                if int(line[22:26]) in residues:
                    for new_line in new_lines:
                        if (
                            new_line[22:26] == line[22:26]
                            and new_line[12:16].strip() == line[12:16].strip()
                        ):
                            line = new_line
            new_pdb.append(line)

        # Write the new PDB file using the basename of the target PDB file
        with open(pdb_tar_name + "_mod.pdb", "w") as f:
            for line in new_pdb:
                f.write(line)
        repeat = input("Would you like to modify another PDB file? (y/n): ")
        if repeat == "n":
            break


if __name__ == "__main__":
    modify_res_pos()
