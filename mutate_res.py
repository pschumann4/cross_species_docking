"""
 Purpose: This script will perform point mutations in a PDB file by replacing specified           
 residues with ALA.                                                                               
 This function will output a new PDB file with the mutated residues.                              
"""

import os

def mutate_res():
    print("\nThis function will perform point mutations in a "
          "PDB file by replacing specified residues with ALA.\n")
    # Ask the user for the working directory
    wd = input('Enter the path to the directory containing the '
               'PDB files you would like to mutate: ')
    # check that the directory exists
    while not os.path.isdir(wd):
        wd = input('Directory not found. Please enter the path to the '
                   'directory containing the PDB files you would like to mutate: ')
    # Change the working directory
    os.chdir(wd)
    # Ask the user for a comma-separated list of residues to mutate
    residues = input('Enter a comma-separated list of residues to mutate: ')
    # Convert to a list by splitting at commas and removing spaces
    residues = residues.split(',')
    for i in range(len(residues)):
        residues[i] = residues[i].strip()
    # Loop through the PDB files in the working directory
    for pdb in os.listdir():
        # Read the lines of the PDB file
        with open(pdb, 'r') as f:
            lines = f.readlines()
            # For each residue number in the list of residues, remove all atoms 
            # that are not either N, CA, C, O, or CB
            for res in residues:
                for i in range(len(lines)):
                    if lines[i][22:26].strip() == res:
                        lines[i] = lines[i][:17] + 'ALA' + lines[i][20:]
                        if lines[i][12:16].strip() not in ['N', 'CA', 'C', 'O', 'CB']:
                            lines[i] = ''
            # Write the new PDB file using the same name as the original 
            # but with '_mutated' appended
            with open(pdb[:-4] + '_mutated.pdb', 'w') as f:
                f.writelines(lines)
    # Move all of the mutated structures to a new folder called 'mutated'
    try:
        os.mkdir('mutated')
    except FileExistsError:
        pass
    for pdb in os.listdir():
        if '_mutated' in pdb:
            # Move the file to the new folder
            os.rename(pdb, 'mutated/' + pdb)
    print('\nDone! All mutated structures have been moved to a '
          'new folder called "mutated".')

if __name__ == '__main__':
    mutate_res()
