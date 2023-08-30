"""
 This script will read a PDB file and return the centroid coordinates from a user-supplied list   
 of residues. The output will be written to a text file with the X, Y, and Z coordinates of the   
 centroid for each binding site.                                                                  
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


def get_centroid(pdb_file):
    """
    Calculates the centroid of a set of coordinates representing a binding site.
    The user will be prompted to enter a comma-separated list of residues.
    """
    # Ask the user for a comma-separated list of residues
    residues = input(
        "Enter a comma-separated list of binding site "
        'residue numbers (e.g., "343, 345, 367" or if chain '
        'info is present "A_343, A_345, A_367"): '
    )
    # Convert the comma-separated list of residues to a list of residues
    if ", " in residues:
        residues = residues.split(", ")
    elif "," in residues:
        residues = residues.split(",")
    if "_" in residues[0]:
        residues = [int(residue.split("_")[1]) for residue in residues]
    else:
        pass
    # Read the lines in the PDB file and extract the coordinates of the
    # specified residues
    coords = []
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM") and int(line[22:26]) in residues:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    # Convert the list to a NumPy array
    coords = np.array(coords)
    # Calculate the centroid of the coordinates
    centroid = np.mean(coords, axis=0)
    # Return the centroid coordinates
    return centroid


# Create a function that will read a PDB file and return the centroid
# coordinates from a user-supplied list of residues
def get_bindingsite_coords():
    """
    This is the main function that will be called by the user.
    """
    # Ask the user for the working directory
    pwd = input("Enter the directory containing the PDB models: ")
    if check_exit(pwd):
        return
    # Check that the working directory exists
    while not os.path.exists(pwd):
        print("The working directory does not exist.")
        pwd = input('Enter the directory (or type "exit"): ')
        if check_exit(pwd):
            return
    os.chdir(pwd)
    # Ask the user for the name of the PDB file
    pdb_file = input("Enter the name of the PDB file: ")
    if pdb_file.endswith(".pdb"):
        pass
    else:
        pdb_file += ".pdb"
    # Check that the PDB file exists
    while not os.path.exists(pdb_file):
        print("The PDB file does not exist.")
        pdb_file = input('Enter the name of the PDB file (or type "exit"): ')
        if pdb_file.endswith(".pdb"):
            pass
        else:
            pdb_file += ".pdb"
        if check_exit(pdb_file):
            return
    # Ask how many binding sites the user wants to analyze
    num_sites = input("How many binding sites do you want to analyze? ")
    if check_exit(num_sites):
        return
    # Use the get_centroid() function to get the centroid coordinates
    # for each binding site
    site_coords = {}
    for i in range(int(num_sites)):
        centroid = get_centroid(pdb_file)
        print("The centroid coordinates for binding site {} are:".format(i + 1))
        print("X: {:.3f}".format(centroid[0]))
        print("Y: {:.3f}".format(centroid[1]))
        print("Z: {:.3f}".format(centroid[2]))
        print()
        site_coords[i + 1] = centroid
    # Write the centroid coordinates to a text file
    pdb_basename = os.path.basename(pdb_file)
    pdb_basename = pdb_basename.split(".")[0]
    with open("{}_centroid_coords.txt".format(pdb_basename), "w") as f:
        for key, value in site_coords.items():
            f.write("Binding site {}:\n".format(key))
            f.write("X: {:.3f}\n".format(value[0]))
            f.write("Y: {:.3f}\n".format(value[1]))
            f.write("Z: {:.3f}\n".format(value[2]))
            f.write("\n")
    print("The centroid coordinates have been written to 'centroid_coords.txt'.")


if __name__ == "__main__":
    get_bindingsite_coords()
