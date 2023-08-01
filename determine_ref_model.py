######################################################################################
# This function will determine the best reference model for a set of PDB files.      #
# It does so by calculating the pairwise RMSD values between the ligand coordinates  #
# in each PDB file. The PDB file with the lowest average RMSD value is suggested.    #
######################################################################################

############################
# Import required packages #
############################

import os
import numpy as np

######################
# Load all functions #
######################

def check_exit(input_str):
    if isinstance(input_str, str):
        input_str = input_str.lower()
        if input_str == 'exit':
            return True
        return False

def determine_ref_model():
    # Ask the user for the working directory
    pwd = input('Enter the directory containing the PDB models: ')
    if check_exit(pwd):
        return
    while not os.path.exists(pwd):
        print('The working directory does not exist.')
        pwd = input('Enter the directory (or type "exit"): ')
        if check_exit(pwd):
            return
    os.chdir(pwd)
    # Ask the user for the ligand name
    ligand = input('Enter the ligand ID as it is found within the PDB files: ')
    if check_exit(ligand):
        return
    # Loop through all of the PDB files in the working directory and check that they contain the ligand
    pdb_files = []
    no_lig = []
    for file in os.listdir():
        if file.endswith('.pdb') and 'model' in file:
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith('HETATM') and line[17:20].strip() == ligand:
                        pdb_files.append(file)
                        break
                    if line.startswith('HETATM') and line[17:20].strip() != ligand:
                        if file not in no_lig:
                            no_lig.append(os.path.basename(file))
    # Print the names of the PDB files that do not contain the ligand
    if len(no_lig) > 0:
        print('WARNING: The following PDB files do not contain the ligand: ')
        for file in no_lig:
            print(file)
    ligand_coords = {}
    # Iterate through the PDB files to get ligand coordinates
    for file in pdb_files:
        coords = {}
        with open(file, 'r') as f:
            for line in f:
                if line.startswith('HETATM') and line[17:20].strip() == ligand:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atom = line[12:16].strip()
                    coords[atom] = [x, y, z]
        ligand_coords[file] = coords
    # Calculate all of the pairwise RMSD values
    rmsd = {}
    for file1 in pdb_files:
        rmsd[file1] = {}
        for file2 in pdb_files:
            if file1 == file2:
                continue
            else:
                # Convert the dictionary of coordinates to a list of coordinates
                coords1 = []
                coords2 = []
                for atom in ligand_coords[file1]:
                    coords1.append(ligand_coords[file1][atom])
                for atom in ligand_coords[file2]:
                    coords2.append(ligand_coords[file2][atom])
                # Calculate the RMSD value
                diff = np.array(coords1) - np.array(coords2)
                rmsd[file1][file2] = np.sqrt(np.sum(diff**2) / len(coords1))
    # Calculate the average RMSD value
    avg_rmsd = {}
    for file in pdb_files:
        avg_rmsd[file] = np.mean(list(rmsd[file].values()))
    # Determine the PDB file with the lowest average RMSD value
    ref_model = min(avg_rmsd, key=avg_rmsd.get)
    print('The suggested reference model is: ' + ref_model)

determine_ref_model()