"""
 This function loops the prep_receptor4.py and prepare_flexreceptor.py 
 scripts over a specified directory of PDB files to create the .pdbqt files 
 for the rigid and flexible residues (if flexible docking is desired).                                  
 requires: prep_receptor4.py, prepare_flexreceptor.py, flex_residues.txt, 
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
    pdb_files = [i for i in pdb_files if i.endswith(".pdb")]

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

        if not os.path.exists("prepare_flexreceptor.py"):
            flex_script = input(
                "Enter the path to the prepare_flexreceptor.py script: "
            )
            flex_script = flex_script.replace('"', "")
            while not os.path.exists(flex_script):
                flex_script = input(
                    "That path does not appear to exist.\nPlease enter the path to "
                    "the prepare_flexreceptor.py script: "
                )
        else:
            flex_script = "prepare_flexreceptor.py"

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
        if not flex_script.startswith('"') and not flex_script.endswith('"'):
            flex_script = '"' + flex_script + '"'
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

        # Prepare the flexible residues for each PDBQT file
        for pdbqt in pdbqt_files:
            # temp_pdbqt is a list of the lines in the PDBQT file that are
            # within +/- 5 residues of the flexible residues
            temp_pdbqt = []
            # res_temp is a list of the res range to subset the PDBQT file
            res_temp = []
            pdbqt_name = pdbqt.split(".pdbqt")[0]
            print("\nPreparing " + pdbqt_name + "...")
            for line in res_list:
                if line.split(":")[0] == pdbqt_name:
                    res = line.split("[")[1].split("]")[0]
                    # If res is empty, skip the file
                    if res == "":
                        print("No flexible residues for " + pdbqt_name)
                        no_flex.append(pdbqt_name)
                        continue
                    # Subset the PDBQT to include +/- 5 residues based on the
                    # first and last residue in the flexible residues list
                    res_id = res.split(",")
                    res_id = [i[3:] for i in res_id] # Remove the res name
                    res_id = sorted(res_id, key=lambda x: int(x))
                    res_id.insert(0, str(int(res_id[0]) - 5))
                    res_id.append(str(int(res_id[-1]) + 5))
                    for i in range(int(res_id[0]), int(res_id[-1]) + 1):
                        res_temp.append(str(i))
                    with open(pdbqt, "r") as f:
                        pdbqt_lines = f.readlines()
                        for line in pdbqt_lines:
                            if line[22:26].strip() in res_temp:
                                temp_pdbqt.append(line)
                    # The subsetted PDBQT file is written to a temp file so that it
                    # can be used as the input for the prepare_flexreceptor.py script
                    with open(pdbqt_name + "_temp.pdbqt", "w") as f:
                        for line in temp_pdbqt:
                            f.write(line)
                    res = res.replace(",", "_")
                    cmd = "python {} -r {}_temp.pdbqt -s {} -g {}_rigid_residues.pdbqt -x {}_flex_residues.pdbqt -v".format(flex_script, pdbqt_name, res, pdbqt_name, pdbqt_name)
                    os.system(cmd)
                    if (
                        not os.path.exists(pdbqt_name + "_flex_residues.pdbqt")
                        and not os.path.exists(pdbqt_name + "_rigid_residues.pdbqt")
                        and pdbqt_name not in no_flex
                    ):
                        failed.append(pdbqt_name)
                        shutil.move(pdbqt, "failed_preps\\" + pdbqt)

        # Remove all of the temp_pdbqt files
        temp_files = os.listdir()
        temp_files = [f for f in temp_files if f.endswith("_temp.pdbqt")]
        for f in temp_files:
            os.remove(f)
        
        # Create a new rigid receptor file by merging the original PDBQT file
        # with the rigid residues file
        pdbqt_files = os.listdir()
        for pdbqt in pdbqt_files:
            new_file = []
            if pdbqt.endswith("modified.pdbqt"):
                print("\nUpdating rigid receptor file for " + pdbqt.split("_", 1)[0] + "...")
                with open(pdbqt, "r") as f:
                    modified = f.readlines()
                # Match the name of the PDBQT file to the name of the rigid residues file
                species = pdbqt.split("_", 1)[0]
                rigid_file = [f for f in pdbqt_files if f.startswith(species) and f.endswith("rigid_residues.pdbqt")]
                for file in rigid_file:
                    with open(file, "r") as f:
                        rigid_pqbqt = f.readlines()
                        rigid_residues = [i[22:26].strip() for i in rigid_pqbqt if i.startswith("ATOM")]
                    # Remove the lines of the residues of the modified PDBQT file that are in the rigid residues file
                    modified = [i for i in modified if i[22:26].strip() not in rigid_residues]
                    # Add all the remaining lines from the modified PDBQT file and the rigid residues file to a new list
                    for line in modified:
                        new_file.append(line)
                    with open(rigid_file[0], "r") as f:
                        rigid = f.readlines()
                        for line in rigid:
                            new_file.append(line)
                    # Overwrite the original rigid residues file with the new file
                    with open(rigid_file[0], "w") as f:
                        for line in new_file:
                            f.write(line)

        # Print the list of failed preps
        if len(failed) > 0:
            print("\nThe following PDBQT files failed to prep:")
            for fail in failed:
                print(fail)
            print("\nConsider a manual prep using the AutoDockTools GUI.")
        else:
            print("\nAll PDBQT files were successfully prepped.")
    else:
        pass


if __name__ == "__main__":
    prep_receptors()
