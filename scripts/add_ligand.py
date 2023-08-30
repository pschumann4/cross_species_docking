"""
 This script is used to generate PDB formatted models from the PDBQT files generated by AutoDock  
 Vina. The script will add the ligand to the PDB file and correct the atom numbers and chain      
 identifiers. If there were flexible residues in the original PDB file, the script will add the   
 flexible residues. If there are multiple models, the script will add the ligand to each model.   
 The script will also move all of the models to a new folder called 'models' in the working       
 directory.                                                                                       
"""

import os
import random
import shutil


def add_flex_residues(file_directory):
    """
    Adds the flexible residues from the Vina output to the original PDB file.
    """
    # List PDB files in directory
    pdb_files = [file for file in os.listdir(file_directory) if file.endswith(".pdb")]

    for pdb_file in pdb_files:
        # Extract PDB file name from path and use it to create output prefix
        pdb_name = os.path.splitext(pdb_file)[0]
        pdb_name = pdb_name.split("_")[0]
        output_prefix = pdb_name + "_flex"
        # Create empty dictionary to store new file lines
        new_lines = {}
        # Start a counter for the output file name
        counter = 1
        # Get list of files in directory
        files = [i for i in os.listdir(file_directory) if i.endswith(".pdbqt") and "flex" in i and pdb_name in i]

        # Loop through the flex PDBQT files in directory that match the PDB file name
        for file in files:
            # Get file path
            file_path = os.path.join(file_directory, file)
            # Open new file and store lines in dictionary using residue number
            # and atom name as key
            with open(file_path, "r") as f:
                for line in f:
                    if line.startswith("ATOM"):
                        residue_name = line[17:20].strip()
                        residue_num = line[22:26].strip()
                        atom_name = line[12:16].strip()
                        new_lines[(residue_name, residue_num, atom_name)] = line
            # Construct output file path
            output_file = os.path.join(
                file_directory, output_prefix + str(counter) + ".pdb"
            )

            # Open original file and write its lines and new lines to output file
            with open(os.path.join(file_directory, pdb_file), "r") as f:
                with open(output_file, "w") as o:
                    for line in f:
                        if line.startswith("ATOM"):
                            residue_name = line[17:20].strip()
                            residue_num = line[22:26].strip()
                            atom_name = line[12:16].strip()
                            if (residue_name, residue_num, atom_name) in new_lines:
                                # Replace the coordinates of the original
                                # line with the coordinates of the new line
                                new_line = list(line)
                                new_line[30:38] = new_lines[
                                    (residue_name, residue_num, atom_name)
                                ][30:38]
                                new_line[38:46] = new_lines[
                                    (residue_name, residue_num, atom_name)
                                ][38:46]
                                new_line[46:54] = new_lines[
                                    (residue_name, residue_num, atom_name)
                                ][46:54]
                                o.write("".join(new_line))
                            else:
                                o.write(line)
                        else:
                            o.write(line)
            # Increment counter
            counter += 1


def add_ligand_to_flex(file_directory):
    """
    Adds the Vina ligand pose to the flexible PDB files.
    """
    # List files in directory
    files = os.listdir(file_directory)

    # Get PDB file info from directory
    for pdb_file in files:
        if pdb_file.endswith(".pdb") and "model" not in pdb_file:
            # Extract PDB file name and number from path
            pdb_name = os.path.splitext(pdb_file)[0]
            pdb_number = pdb_name[-1]
            pdb_name = pdb_name.split("_")[0]

            # Loop through files in directory that contain the word "ligand"
            # in the file name and match the PDB file name and number
            for file in files:
                if (
                    file.endswith(pdb_number + ".pdbqt")
                    and "ligand" in file
                    and pdb_name in file
                ):
                    # Get file path
                    file_path = os.path.join(file_directory, file)
                    # Create output file path
                    output_file = os.path.join(
                        file_directory, pdb_name + "_flex_model" + pdb_number + ".pdb"
                    )
                    # Create empty list to store lines
                    lines = []

                    # Open new file and store the 'HETATM' lines in list
                    with open(file_path, "r") as f:
                        for line in f:
                            if line.startswith("HETATM"):
                                # Remove the segment identifier [70:76] from the line
                                # and add 6 space characters in its place
                                line = line[:70] + " " * 6 + line[76:]
                                # Remove the active or inactive flags from the
                                # element symbol [78:80]
                                line = line[:78] + line[79:]
                            # Store the line in the list unless it's a hydrogen
                            if line[76:78].strip() != "H":
                                lines.append(line)

                    # Open PDB file and store the lines in list
                    with open(os.path.join(file_directory, pdb_file), "r") as f:
                        for line in f:
                            lines.append(line)

                    # Open output file and write the lines to it
                    with open(output_file, "w") as o:
                        for line in lines:
                            o.write(line)


def add_ligand_to_rigid(file_directory):
    """
    Adds the Vina ligand pose to the rigid PDB files.
    """
    # List files in directory
    files = os.listdir(file_directory)
    # List file names to track which files are flexible
    flex_files = []

    # Get the flexible PDB files in directory
    for file in files:
        if file.endswith(".pdb") and "flex" in file:
            # Extract PDB file name from path
            pdb_name = os.path.splitext(file)[0]
            pdb_name = pdb_name.split("_")[0]
            flex_files.append(pdb_name)

    # Get all PDB files in directory
    for file in files:
        if file.endswith(".pdb"):
            # Extract PDB file name from path
            pdb_name = os.path.splitext(file)[0]
            pdb_name = pdb_name.split("_")[0]

            # Check if the file name is in the list of flexible files
            # and if it is, skip it
            if pdb_name not in flex_files and "model" not in file:
                # Loop through files in directory that contain the word
                # "ligand" in the file name and match the PDB file name
                for pdb_file in files:
                    if "ligand" in pdb_file and pdb_name in pdb_file:
                        # Get file path
                        file_path = os.path.join(file_directory, pdb_file)
                        # Get the file number
                        file_number = os.path.splitext(pdb_file)[0].split("_")[-1]
                        # Create output file path
                        output_file = os.path.join(
                            file_directory,
                            pdb_name + "_rigid_model" + file_number + ".pdb",
                        )
                        # Create empty list to store lines
                        lines = []

                        # Open new file and store the 'HETATM' lines in list
                        with open(file_path, "r") as f:
                            for line in f:
                                if line.startswith("HETATM"):
                                    # Remove the segment identifier from the line
                                    line = line[:70] + " " * 6 + line[76:]
                                    # Remove the active or inactive flags from
                                    # the element symbol
                                    line = line[:78] + line[79:]
                                # Store the line in the list unless it's a hydrogen
                                if line[76:78].strip() != "H":
                                    lines.append(line)

                        # Open PDB file and store the lines in list
                        with open(os.path.join(file_directory, file), "r") as f:
                            for line in f:
                                lines.append(line)

                        # Open output file and write the lines to it
                        with open(output_file, "w") as o:
                            for line in lines:
                                o.write(line)
            else:
                continue


def check_chain_ids(pdb_file, ref_pdb, ligand):
    """
    Checks if the chain identifiers in the PDB file
    match the chain identifiers in the reference PDB file.
    """
    # Initialize list to store the modified records
    modified_records = []

    # Open the ref_pdb and extract the chain identifier for ATOMS and HETATMS
    ref_chain_id = ""
    lig_chain_id = ""
    with open(ref_pdb, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                ref_chain_id = line[21:22]
            if line.startswith("HETATM") and line[17:20].strip() == ligand:
                lig_chain_id = line[21:22]

    # Open the pdb_file
    with open(pdb_file, "r") as f:
        # Read each line in the file
        for line in f:
            # Check if the line is an ATOM record
            if line.startswith("ATOM"):
                # Check if the chain identifier is the same as the reference
                # chain identifier
                if line[21:22] != ref_chain_id:
                    modified_record = line[:21] + ref_chain_id + line[22:]
                    modified_records.append(modified_record)
                else:
                    modified_records.append(line)

            # Check if the line is a HETATM record
            if line.startswith("HETATM"):
                # Check if the chain identifier is the same as the reference
                # ligand chain identifier
                if line[21:22] != lig_chain_id:
                    modified_record = line[:21] + lig_chain_id + line[22:]
                    modified_records.append(modified_record)
                else:
                    modified_records.append(line)

            # If the line is not an ATOM or HETATM record, add it to the
            # list of modified records
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                modified_records.append(line)

    # Write the modified records to the same file
    with open(pdb_file, "w") as f:
        for record in modified_records:
            f.write(record)


def move_hetatms(pdb_file):
    """
    Reorganizes the records in the PDB file so that the HETATM records 
    are at the beginning of the file. This is necessary
    for the PDB files to be read correctly by the 'PLIP' program.
    """
    # Initialize lists to store records
    hetatms = []
    atoms = []
    end_lines = []
    other_lines = []

    # Open the PDB file and read the lines
    with open(pdb_file, "r") as f:
        for line in f:
            # Check if the line starts with 'HETATM'
            if line.startswith("HETATM"):
                # Check that the element symbol matches the atom name
                if line[13:14].strip() != line[76:78].strip():
                    # Replace the element symbol with the atom name
                    line = line[:77] + line[13:14].strip() + line[78:]
                # Add the line to the 'HETATM' list
                hetatms.append(line)
                # Reorder the hetatm records based on the atom name [12:16]
                # in alphabetical order
                hetatms.sort(key=lambda x: x[12:16])
            # Check if the line starts with 'ATOM'
            elif line.startswith("ATOM"):
                # Add the line to the 'ATOM' list
                atoms.append(line)
            # Check if the line starts with 'TER' or 'END'
            elif line.startswith("TER") or line.startswith("END"):
                # Add the line to the 'end_lines' list
                end_lines.append(line)
            # Check for 'CRYST1' and 'REMARK' records
            elif line.startswith("CRYST1") or line.startswith("REMARK"):
                # Add the line to the 'other_lines' list
                other_lines.append(line)
    # Combine all of the lists into a single list
    all_lines = other_lines + hetatms + atoms + end_lines
    # Correct the atom number order
    # Create a counter to keep track of the atom number
    counter = 1

    for i, line in enumerate(all_lines):
        # Check if the line starts with 'HETATM' or 'ATOM'
        if line.startswith("HETATM") or line.startswith("ATOM"):
            # Extract the atom number from the line
            atom_number = int(line[6:11])
            # Check if the atom number is not equal to the current index
            if atom_number != counter:
                # Replace the atom number with the current index
                all_lines[i] = line[:6] + str(counter).rjust(5) + line[11:]
            # Increment the counter
            counter += 1
        # Remove all the 'CONECT' and 'ROOT' and 'ENDROOT' records
        if (
            line.startswith("CONECT")
            or line.startswith("ROOT")
            or line.startswith("ENDROOT")
        ):
            all_lines[i] = ""

    # Open the PDB file and write the lines back to the file
    with open(pdb_file, "w") as f:
        for line in all_lines:
            f.write(line)


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


def add_ligand():
    """
    Main function
    """
    # Prompt user for directory path
    pwd = input(
        "Enter the directory containing the PDB files and related PDBQT files: "
    )

    # Check that the working directory exists
    while not os.path.exists(pwd):
        print("The working directory does not exist.")
        pwd = input('Enter the directory (or type "exit"): ')
        if check_exit(pwd):
            return

    # Set the current working directory to the directory path
    os.chdir(pwd)
    
    ref = ""
    ref_name = ""

    # Search the directory for a file that starts with "ref_"
    for file in os.listdir(pwd):
        if file.startswith("ref_"):
            ref = input("Is {} the reference PDB file? (y/n): ".format(file))
            ref = ref.lower()
            # check for valid input
            while ref.lower() not in ["y", "n"]:
                ref = input("Please enter y or n: ")
            if ref == "y":
                ref_name = file
                break
    if ref == "n":
        ref_name = input("Enter the name of the reference PDB file (if none, type 'random'): ")
        if check_exit(ref_name):
            return
        if ref_name.endswith(".pdb") or ref_name == "random":
            pass
        else:
            ref_name += ".pdb"
        if ref_name == "random":
            # Create a list of PDB files in the directory
            pdb_files = [f for f in os.listdir(pwd) if f.endswith(".pdb")]
            # Randomly select a PDB file from the list
            ref_name = random.choice(pdb_files)
            print("Randomly selected reference PDB file: {}".format(ref_name))
    
    # Check that the reference file exists and if not ask for it again
    while not os.path.exists(ref_name):
        print("The reference file does not exist.")
        ref_name = input("Enter the name of the reference PDB file: ")
        if check_exit(ref_name):
            return
        if ref_name.endswith(".pdb"):
            pass
        else:
            ref_name += ".pdb"

    # Ask for the name of the ligand
    ligand = input("Enter the ligand ID as it is found within the reference PDB: ")
    if check_exit(ligand):
        return

    # Inform the user that flexible residues are being added to the PDB files
    print("Adding flexible residues to PDB files, if any...")
    add_flex_residues(pwd)
    # Create a list of file names to track which files are flexible
    flex_files = []

    # Loop through the files in the directory
    for file in os.listdir(pwd):
        # Check if the file is a PDB file and if it is a flexible model
        if file.endswith(".pdb") and "flex" in file:
            # Get the base name of the file and append it to
            # the list of flexible files
            file_name = os.path.splitext(file)[0]
            file_name = file_name.split("_")[0]
            flex_files.append(file_name)

    # Re-loop through the files in the directory
    for file in os.listdir(pwd):
        # Check if the file is a PDB file and if it is a flexible model
        if file.endswith(".pdb"):
            # Get the base name of the file
            file_name = os.path.splitext(file)[0]
            file_name = file_name.split("_")[0]
            # Check if the file is in the list of flexible files
            if file_name in flex_files and "flex" in file:
                add_ligand_to_flex(pwd)
                print("Adding ligand to flexible model: " + file)
            # If not in flexible files or the reference model,
            # then add ligand to rigid model
            elif (
                file_name not in flex_files and "flex" not in file and file != ref_name
            ):
                add_ligand_to_rigid(pwd)
                print("Adding ligand to rigid model: " + file)

    # Re-loop through the files in the directory to apply the move_hetatms()
    # and check_chain_identifiers() functions
    for file in os.listdir(pwd):
        if file.endswith(".pdb") and "model" in file:
            print("Correcting for ligand addition to PDB file " + file)
            move_hetatms(file)
            check_chain_ids(file, ref_name, ligand)
        # Remove all PDB files that contain that word 'flex'
        # but do not contain the word 'model'
        if file.endswith(".pdb") and "flex" in file and "model" not in file:
            os.remove(file)
    print("Correcting reference PDB file...")
    move_hetatms(ref_name)

    # Move all of the models to a new folder called 'models'
    print(
        "All PDB models are being moved to the 'models' "
        "folder in the working directory..."
    )
    if not os.path.exists("models"):
        os.mkdir("models")
    for file in os.listdir(pwd):
        if file.endswith(".pdb") and "model" in file:
            shutil.move(file, "models")

    # Move a copy of the reference PDB file to the 'models' folder
    shutil.copy(ref_name, "models")
    # Check that the 'models' folder exists and is not empty
    if os.path.exists("models") and os.listdir("models"):
        print("DONE!")
    else:
        print('ERROR: Models could not be moved to the "models" folder.')

if __name__ == "__main__":
    add_ligand()