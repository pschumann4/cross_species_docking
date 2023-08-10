"""
 This function loops the prepare_flexreceptor.py script over a specified 
 directory of PDB files to create the .pdbqt files for the rigid 
 and flexible residues.                                  
 requires: prepare_flexreceptor.py script, flex_residues.txt file, 
 and a directory of PDB files   
 returns: .pdbqt files for the rigid and flexible residues for each PDB 
 file in the directory     
 The prepare_flexreceptor.py script needs to be ran in python v2.7 and the 
 user must have MGL Tools intalled on their computer. "python" is the default 
 command for python v2.7 once MGL Tools is installed.
"""

import os


def prep_flex_receptors():
    """
    Main function
    """
    # Get the path to the directory containing the PDBQT files
    pdbqt_dir = input("Enter the path to the directory containing the PDBQT files: ")
    # Check if the path exists and if not, ask for it again
    while not os.path.exists(pdbqt_dir):
        pdbqt_dir = input(
            "That path does not appear to exist.\nPlease enter the path to "
            "the directory containing the PDBQT files: "
        )
    # Change the working directory to the directory containing the PDBQT files
    os.chdir(pdbqt_dir)
    # Get the list of PDBQT files in the directory
    pdbqt_files = os.listdir(pdbqt_dir)

    # Specify name of prepare_flexreceptor.py script
    script = input("Enter the path to the prepare_flexreceptor.py script: ")
    # Remove the "" from the path if they were included
    if script.startswith('"') and script.endswith('"'):
        script = script[1:-1]
    # Check if the path exists and if not, ask for it again
    while not os.path.exists(script):
        script = input(
            "That path does not appear to exist.\nPlease enter the path to "
            "the prepare_flexreceptor.py script: "
        )
        if script.startswith('"') and script.endswith('"'):
            script = script[1:-1]

    # Specify name of flex_residues.txt file
    residues = input("Enter the path to the flex_residues.txt file: ")
    # Remove the "" from the path if they were included
    if residues.startswith('"') and residues.endswith('"'):
        residues = residues[1:-1]
    # Check if the path exists and if not, ask for it again
    while not os.path.exists(residues):
        residues = input(
            "That path does not appear to exist.\nPlease enter the path to "
            "the flex_residues.txt file: "
        )
        if residues.startswith('"') and residues.endswith('"'):
            residues = residues[1:-1]

    # Open the flex_residues.txt file and get the list of PDBQT files
    # and their flexible residues
    with open(residues, "r") as f:
        res_list = f.readlines()
        # Remove header and footer lines
        res_list = res_list[2:-2]

    # Loop the prepare_flexreceptor.py script over the PDBQT files
    for pdbqt in pdbqt_files:
        if pdbqt.endswith(".pdbqt"):
            # Get the name of the PDBQT file
            pdbqt_name = pdbqt.split(".pdbqt")[0]
            # Match the PDBQT file name to the PDBQT ID in the flex_residues.txt
            # file and get the residues to set as flexible
            for line in res_list:
                if line.split(":")[0] == pdbqt_name:
                    res = line.split("[")[1].split("]")[0]
                    # Create the command to run the prepare_flexreceptor.py script
                    cmd = (
                        "python "
                        + script
                        + " -r "
                        + pdbqt_name
                        + ".pdbqt -s "
                        + res
                        + " -g "
                        + pdbqt_name
                        + "_rigid_residues.pdbqt -x "
                        + pdbqt_name
                        + "_flex_residues.pdbqt -v"
                    )
                    os.system(cmd)
        else:
            continue
    print("DONE!")


if __name__ == "__main__":
    prep_flex_receptors()
