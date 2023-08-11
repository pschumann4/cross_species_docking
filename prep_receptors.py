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
import shutil


def prep_receptors():
    """
    Main function
    """
    # Get the path to the directory containing the PDBQT files
    pdb_dir = input("Enter the path to the directory containing the PDB files: ")
    # Check if the path exists and if not, ask for it again
    while not os.path.exists(pdb_dir):
        pdb_dir = input(
            "That path does not appear to exist.\nPlease enter the path to "
            "the directory containing the PDB files: "
        )
    # Change the working directory to the directory containing the PDBQT files
    os.chdir(pdb_dir)
    # Get the list of PDB files in the directory
    pdb_files = os.listdir(pdb_dir)
    pdb_files = [f for f in pdb_files if f.endswith(".pdb")]

    # Prepare PDBQT files?
    prep_rec = input("\nDo you want to prepare the PDBQT files? (y/n): ").lower()
    while prep_rec != "y" and prep_rec != "n":
        prep_rec = input("Please enter y or n: ").lower()

    if prep_rec == "y":
        # Specify name of prepare_receptor4.py script
        if not os.path.exists("prep_receptor4.py"):
            prep_pdbqt = input("Enter the path to the prep_receptor4.py script: ")
            prep_pdbqt = prep_pdbqt.replace('"', "")
            while not os.path.exists(prep_pdbqt):
                prep_pdbqt = input(
                    "That path does not appear to exist.\nPlease enter the path to "
                    "the prep_receptor4.py script: "
                )
        else:
            prep_pdbqt = "prep_receptor4.py"
        if not prep_pdbqt.startswith('"') and not prep_pdbqt.endswith('"'):
            prep_pdbqt = '"' + prep_pdbqt + '"'
        failed_preps = []
        # Create a folder called "pdbqt_files" if it does not exist
        if not os.path.exists("pdbqt_files"):
            os.mkdir("pdbqt_files")
        # Run the prepare_receptor4.py script on each PDB file in the directory
        for pdb in pdb_files:
            cmd = "python {} -r {} -v -o {} -A checkhydrogens".format(prep_pdbqt, pdb, pdb[:-4] + ".pdbqt")
            os.system(cmd)
            shutil.move(pdb[:-4] + ".pdbqt", "pdbqt_files\\" + pdb[:-4] + ".pdbqt")
        # Check for failed preps
        pdbqt_files = os.listdir("pdbqt_files")
        for pdb in pdb_files:
            if pdb[:-4] + ".pdbqt" not in pdbqt_files:
                failed_preps.append(pdb)
        if len(failed_preps) > 0:
            print("\nThe following files were not converted to PDBQT:")
            for pdb in failed_preps:
                print(pdb)
            print("\nConsider a manual prep using the MGLTools GUI.")
        else:
            print("\nAll PDB files were converted to PDBQT.")
    else:
        pass

    # Prep flexible residues?
    prep_flex = input("\nDo you want to prepare the flexible residues? (y/n): ").lower()
    while prep_flex not in ["y", "n"]:
        prep_flex = input("Please enter y or n: ").lower()

    # If the user wants to prepare the flexible residues
    if prep_flex == "y":
        # Change the working directory to the directory containing the PDBQT files
        os.chdir("pdbqt_files")
        pdbqt_files = os.listdir()
        pdbqt_files = [f for f in pdbqt_files if f.endswith(".pdbqt")]

        if not os.path.exists("prepare_flexreceptors.py"):
            flex_script = input(
                "Enter the path to the prepare_flexreceptors.py script: "
            )
            flex_script = flex_script.replace('"', "")
            while not os.path.exists(flex_script):
                flex_script = input(
                    "That path does not appear to exist.\nPlease enter the path to "
                    "the prepare_flexreceptors.py script: "
                )
        else:
            flex_script = "prepare_flexreceptors.py"

        if not flex_script.startswith('"') and not flex_script.endswith('"'):
            flex_script = '"' + flex_script + '"'

        # Check if the file "flex_residues.txt" exists in the directory
        if not os.path.exists("flex_residues.txt"):
            residues = input("Enter the path to the flex_residues.txt file: ")
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
        else:
            residues = "flex_residues.txt"
        # Read the flex_residues.txt file
        with open(residues, "r") as f:
            res_list = f.readlines()
            # Remove header and footer lines
            res_list = res_list[2:-2]
        # Initialize a list to track failed preps
        failed = []
        # Initialize list to track PDBQT files that have no flexible residues
        no_flex = []
        # Create a folder for the failed preps
        if not os.path.exists("failed_preps"):
            os.mkdir("failed_preps")
        # Loop the prepare_flexreceptor.py script over the PDBQT files
        for pdbqt in pdbqt_files:
            # Get the name of the PDBQT file
            pdbqt_name = pdbqt.split(".pdbqt")[0]
            print("\nPreparing " + pdbqt_name + "...")
            # Match the PDBQT file name to the PDBQT ID in the flex_residues.txt
            # file and get the residues to set as flexible
            for line in res_list:
                if line.split(":")[0] == pdbqt_name:
                    res = line.split("[")[1].split("]")[0]
                    # If res is empty, skip the file
                    if res == "":
                        print("No flexible residues for " + pdbqt_name)
                        no_flex.append(pdbqt_name)
                        continue
                    # Replace the commas with underscores
                    res = res.replace(",", "_")
                    # Create the command to run the prepare_flexreceptor.py script
                    cmd = "python {} -r {}.pdbqt -s {} -g {}_rigid_residues.pdbqt -x {}_flex_residues.pdbqt -v".format(flex_script, pdbqt_name, res, pdbqt_name, pdbqt_name)
                    os.system(cmd)
                    # If there are no flex_residues and rigid_residues .pdbqt files, add the PDBQT file name to the failed list
                    if (
                        not os.path.exists(pdbqt_name + "_flex_residues.pdbqt")
                        and not os.path.exists(pdbqt_name + "_rigid_residues.pdbqt")
                        and pdbqt_name not in no_flex
                    ):
                        failed.append(pdbqt_name)
                        # Move the PDB file to the failed_preps folder
                        shutil.move(pdbqt, "failed_preps\\" + pdbqt)
        # Print the list of failed preps
        if len(failed) > 0:
            print("The following PDBQT files failed to prep:")
            for fail in failed:
                print(fail)
            print("Consider a manual prep using the AutoDockTools GUI.")
        else:
            print("All PDBQT files were successfully prepped.")
    else:
        pass


if __name__ == "__main__":
    prep_receptors()
