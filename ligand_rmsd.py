####################################################################################################
# This script will calculate the RMSD between the ligand in a reference PDB file and the ligand in #
# a set of other PDB files.                                                                        #
# The output will be a text file containing the RMSD values for each PDB file.                     #
####################################################################################################

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

def ligand_rmsd():
    # Prompt user for the working directory
    pwd = input('Enter the directory containing the PDB models: ')
    if check_exit(pwd):
        return
    
    # Check that the working directory exists
    while not os.path.exists(pwd):
        print('The working directory does not exist.')
        pwd = input('Enter the directory (or type "exit"): ')
        if check_exit(pwd):
            return
        
    # Change the working directory
    os.chdir(pwd)

    # Prompt user for the reference PDB file
    ref = input('Enter the name of the reference PDB file: ')
    if check_exit(ref):
        return
    
    # Check if the PDB file ends with '.pdb'
    if not ref.endswith('.pdb'):
        ref += '.pdb'

    # Check that the PDB file exists
    while not os.path.exists(ref):
        print('The pdb file you entered does not appear to exist. The entry is case sensitive.')
        ref = input('Enter the name of the reference PDB file (or type "exit"): ')
        if not ref.endswith('.pdb'):
            ref += '.pdb'
        if check_exit(ref):
            return
        
    # Prompt user for the ligand name
    ligand = input('Enter the ligand ID as it is found within the PDB files: ')
    if check_exit(ligand):
        return
    
    # Read the PDB file and store the HETATM lines for the specified ligand in a list
    hetatm1 = []
    with open(ref, 'r') as f:
        for line in f:
            if line.startswith('HETATM') and line[17:20].strip() == ligand:
                hetatm1.append(line)
    # Check that the ligand is present in the PDB file
    if len(hetatm1) == 0:
        print('Error: The ligand was not found in the PDB file!')
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
        if file.endswith('.pdb') and 'model' in file and file != ref:
            # Get the base name of the PDB file
            base = os.path.basename(file)
            base = os.path.splitext(base)[0]
            # Read the PDB file and store the HETATM lines for the specified ligand in a list
            hetatm2 = []
            with open(file, 'r') as f:
                for line in f:
                    if line.startswith('HETATM') and line[17:20].strip() == ligand:
                        hetatm2.append(line)
            # Check that the ligand is present in the PDB file, if not, skip the file
            if len(hetatm2) == 0:
                print('The ligand was not found in ' + file + '.')
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
                    print('ERROR: The atoms in ' + ref + ' and ' + file + ' are not identical.\nThis must be corrected before the RMSD can be calculated properly.')
                    # Print the lists of atoms in each PDB file
                    print('Atoms in ' + ref + ': ' + str(atoms1))
                    print('Atoms in ' + file + ': ' + str(atoms2))
                    return
                # Perform RMSD calculation
                diff = coords1 - coords2
                rmsd = np.sqrt(np.sum(diff ** 2) / len(coords1))
                # Round the rmsd value to 3 decimal places
                rmsd = round(rmsd, 3)
                # Add the RMSD value to the dictionary
                rmsd_dict[base] = rmsd

    # Write the RMSD values to a file
    with open('ligand_rmsd.txt', 'w') as f:
        # Add a title line
        f.write('Ligand RMSD values\n' + '\n')
        # Add the RMSD values for each PDB file
        for pdb2, rmsd in rmsd_dict.items():
            f.write('%s: %s\n' % (pdb2, rmsd))

    # Remove the quotes from the text file and the spaces between the residue name and commas
    with open('ligand_rmsd.txt', 'r') as f:
        lines = f.readlines()
    with open('ligand_rmsd.txt', 'w') as f:
        for line in lines:
            line = line.replace("'", '')
            line = line.replace(', ', ',')
            f.write(line)

    # Check that the text file was created
    if os.path.exists('ligand_rmsd.txt'):
        print('The ligand RMSD values were written to ligand_rmsd.txt in your working directory.')
    return

ligand_rmsd()