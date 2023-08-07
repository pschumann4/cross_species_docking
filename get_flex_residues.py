"""
 This script is used to identify the flexible residues in the binding 
 site of a protein by comparing the binding site structure to a reference structure. 
 Residues are considered flexible if they contain an atom that is greater 
 than 3 Angstroms away from the corresponding atom in the reference structure.
"""

import os

import numpy as np


def euclidean3d(v1, v2):
    """
    Faster implementation of euclidean distance for the 3D case.
    """
    if not len(v1) == 3 and len(v2) == 3:
        return None
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def min_dist(pdb_file, ligand, min_cutoff=5):
    """
    Determines the residues that are within a 5 Angstrom radius of the ligand.
    """
    # Create a list to hold the coordinates of the ligand atoms
    ligand_coords = []
    # Read PDB file
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        # Add the coordinates of the ligand atoms to the ligand_coords list
        if line.startswith("HETATM"):
            if line[17:20].strip() == ligand:
                ligand_coords.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
    # Create a dictionary to store the residue coordinates
    res_coords = {}
    # On a per residue basis, calculate the residue centroids
    for line in lines:
        if line.startswith("ATOM"):
            # Extract the residue number (resnr) from the line
            resnr = line[22:26].strip()
            # Extract the coordinates from the line
            coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            # If the resnr is not in the res_coords dictionary,
            # add it and set the value to a list containing the coordinates
            if not resnr in res_coords:
                res_coords[resnr] = [coords]
            # If the resnr is in the res_coords dictionary, append the coordinates to the list
            else:
                res_coords[resnr].append(coords)
    # Calculate the minimum distance between the each of the ligand and the residue atoms on a per residue basis
    min_dist = {}
    for resnr, coords in res_coords.items():
        min_dist[resnr] = 100
        for ligand_coord in ligand_coords:
            for coord in coords:
                dist = euclidean3d(ligand_coord, coord)
                if dist < min_dist[resnr]:
                    min_dist[resnr] = dist
    # Save only the residue numbers that have a minimum distance less than 5 Angstroms
    min_dist = {resnr: dist for resnr, dist in min_dist.items() if dist < min_cutoff}
    return min_dist


def residue_dist(pdb1, pdb2, dist_cutoff=3):
    """
    Determines the residues that have an atom that is greater than 3 Angstroms
    away from the corresponding atom in the reference structure.
    """
    # Read the PDB files
    with open(pdb1, "r") as f:
        lines1 = f.readlines()
    with open(pdb2, "r") as f:
        lines2 = f.readlines()
    # Create a dictionary to store the residue distances
    dist_dict = {}
    # Iterate through the lines in the PDB file, create temporary lists
    # to store the coordinates of the residues
    for line in lines1:
        # Create a temporary list to store the coordinates of the residue atoms
        coords1 = []
        if line.startswith("ATOM"):
            # Extract the coordinates, residue name, and residue number
            coords1 = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            resname1 = line[17:20].strip()
            resnr1 = line[22:26].strip()
            resid1 = resname1 + resnr1
            # Iterate through the lines in the PDB file, create temporary lists
            # to store the coordinates of the residues
            for line2 in lines2:
                coords2 = []
                if (
                    line2.startswith("ATOM")
                    and line2[22:26].strip() == line[22:26].strip()
                    and line2[13:16].strip() == line[13:16].strip()
                ):
                    # Extract the residue name and number
                    resname2 = line2[17:20].strip()
                    resnr2 = line2[22:26].strip()
                    resid2 = resname2 + resnr2
                    if resid2 == resid1:
                        # Extract the coordinates
                        coords2 = [
                            float(line2[30:38]),
                            float(line2[38:46]),
                            float(line2[46:54]),
                        ]
                        # Calculate the euclidean distance between the two coordinates
                        if coords1 and coords2:
                            dist = euclidean3d(coords1, coords2)
                        # Check if the resid is already in the dictionary,
                        # if not, add it. If there is already a distance for that residue,
                        # retain the larger distance
                        if resid1 in dist_dict:
                            if dist > dist_dict[resid1]:
                                dist_dict[resid1] = dist
                            else:
                                pass
                        else:
                            dist_dict[resid1] = dist
    # Filter the dictionary to only include residues that are greater than the cutoff distance
    filtered_dist = {k: v for k, v in dist_dict.items() if v > dist_cutoff}
    return filtered_dist


def compare_pdbs(pdb1, pdb2, ligand):
    """
    Uses the min_dist and residue_dist functions to determine the top 10 flexible residues
    (i.e., those that are the most distant from the reference structure).
    """
    # Calculate the minimum distance between the ligand and the residue atoms
    # on a per residue basis in the first pdb file
    close_residues = min_dist(pdb1, ligand)
    close_residues = close_residues.keys()
    # Get the list of flexible residues
    flex_residues = residue_dist(pdb1, pdb2)
    # Filter the list of flexible residues to only include residues
    # that are in the close_residues list
    flex_residues_close = {
        residue: rmsd
        for residue, rmsd in flex_residues.items()
        if residue[3:] in close_residues
    }
    # Rank the flexible residues by dist value and return the top 10
    flex_residues10 = sorted(
        flex_residues_close.items(), key=lambda x: x[1], reverse=True
    )[:10]
    # Return a list of the top 10 flexible residues keys
    return [residue[0] for residue in flex_residues10]


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


def get_flex_residues():
    """
    Main function
    """
    # Ask the used for the name of the working directory
    working_dir = input(
        "Enter the name of the file directory containing the PDB files to analyze: "
    )
    if check_exit(working_dir):
        return
    # Set the path to the working directory
    os.chdir(working_dir)
    # Ask the user for the name of the ligand
    ligand = input("Enter the ligand ID as it is found within the reference PDB: ")
    if check_exit(ligand):
        return
    # Ask the user for the reference (pdb1)
    pdb1 = input("Enter the name of the reference structure (with ligand bound): ")
    if not pdb1.endswith(".pdb"):
        pdb1 += ".pdb"
    # Check that the pdb file exists and if not, ask the user to enter it again
    while not os.path.exists(pdb1):
        print("This pdb file does not exist.")
        pdb1 = input(
            "Enter the name of the reference structure (with the ligand bound): "
        )
        if check_exit(pdb1):
            return
        if not pdb1.endswith(".pdb"):
            pdb1 += ".pdb"
    # Check that the pdb file contains the ligand
    with open(pdb1, "r") as f:
        lines = f.readlines()
    ligand_in_pdb = False
    for line in lines:
        if line.startswith("HETATM") and line[17:20].strip() == ligand:
            ligand_in_pdb = True
            break
    # If the ligand is not in the pdb file, ask the user to enter it again
    while not ligand_in_pdb:
        print("This pdb file does not contain the ligand.")
        pdb1 = input(
            "Enter the name of the reference structure "
            'containing the ligand or type "exit": '
        )
        if check_exit(pdb1):
            return
    # Check that the pdb file exists
    # and if not, ask the user to enter the name of the pdb file again
    while not os.path.exists(pdb1):
        print("This pdb file does not exist.")
        pdb1 = input(
            "Enter the name of the reference structure containing "
            'the ligand or type "exit": '
        )
        if check_exit(pdb1):
            return
    # Get the base name of the binding site structure file
    pdb1_base = os.path.basename(pdb1)
    # Get the name of the binding site structure file as the
    # characters before the first underscore
    pdb1_name = pdb1_base.split("_")[0]
    # Get the pdb files that are not the reference binding site
    pdb2s = []
    for file in os.listdir():
        if file.endswith(".pdb") and file != pdb1 and pdb1_name not in file:
            pdb2s.append(file)
    # Check that there are pdb files in the working directory
    # that are not the reference binding site
    if len(pdb2s) == 0:
        print("Error: There are no other pdb files in the working directory.")
        return

    # Create a dictionary to store the flexible residues
    flex_residues = {}
    # Loop through the pdb files and get the flexible residues
    for pdb2 in pdb2s:
        # State which pdb file is being processed
        print("Processing %s" % pdb2 + "...")
        flex_residues[pdb2] = compare_pdbs(pdb1, pdb2, ligand)
        # Check if there are any flexible residues
        if len(flex_residues[pdb2]) == 0:
            print("No flexible residues were found.")
        else:
            print("The flexible residues are: %s" % flex_residues[pdb2])
    # Combine all of the flexible residues into a single list and remove duplicates
    all_flex_residues = list(
        set([res for pdb2, residues in flex_residues.items() for res in residues])
    )
    print("The combined list of flexible residues is: %s" % all_flex_residues)

    # Write the flexible residues to a file
    with open("flex_residues.txt", "w") as f:
        f.write("FLEXIBLE RESIDUES\n" + "\n")
        for pdb2, residues in flex_residues.items():
            f.write("%s: %s\n" % (pdb2.split(".pdb")[0], residues))
        f.write("\nCombined list: %s\n" % all_flex_residues)
    with open("flex_residues.txt", "r") as f:
        lines = f.readlines()
    with open("flex_residues.txt", "w") as f:
        for line in lines:
            line = line.replace("'", "")
            line = line.replace(", ", ",")
            f.write(line)
    # Check if a folder called 'details' exists and if not,
    # create it and move the text file into it
    if not os.path.exists("details"):
        os.mkdir("details")
    os.rename("flex_residues.txt", "details/flex_residues.txt")
    print(
        "\nA summary of these results have been written to the file "
        '"flex_residues.txt" and can be found in the "details" folder.'
    )
    return


if __name__ == "__main__":
    get_flex_residues()
