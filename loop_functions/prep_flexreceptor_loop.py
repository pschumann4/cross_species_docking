import os

def prep_flex_receptors():
    """This function loops the prepare_flexreceptor.py script over a specified directory of PDB files to create 
    the .pdbqt files for the rigid and flexible residues.
    requires: prepare_flexreceptor.py script, flex_residues.txt file, and a directory of PDB files
    returns: .pdbqt files for the rigid and flexible residues for each PDB file in the directory
    """
    
    # Get the path to the directory containing the PDB files
    pdbqt_dir = os.getcwd()
    pdbqt_files = os.listdir(pdbqt_dir)

    # Specify name of prepare_flexreceptor.py script
    script = "prepare_flexreceptor.py"

    # Specify name of flex_residues.txt file
    residues = 'flex_residues.txt'
    with open(residues, 'r') as f:
        res_list = f.readlines()
        # Remove header and footer lines
        res_list = res_list[2:-2]

    # Loop the prepare_flexreceptor.py script over the PDB files
    for pdbqt in pdbqt_files:
        if pdbqt.endswith('.pdbqt'):    
            # Match the PDB file name to the PDB ID in the flex_residues.txt file and get the residues to set as flexible
            for line in res_list:
                if line.split(':')[0] == pdbqt:
                    res = line.split('[')[1].split(']')[0]
                    res = res.replace(',', '_')
                    pdbqt_name = pdbqt.split('.pdbqt')[0]
                    # Create the command to run the prepare_flexreceptor.py script
                    cmd = 'pythonsh ' + script + ' -r ' + pdbqt_dir + '/' + pdbqt + ' -s ' + res + ' -g ' + pdbqt_name + '_rigid_residues.pdbqt -x ' + pdbqt_name + '_flex_residues.pdbqt' + ' -v'
                    os.system(cmd)
        else:
            continue
    print('Done!')
        
if __name__ == '__main__':
    prep_flex_receptors()