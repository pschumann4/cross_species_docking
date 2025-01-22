import os

def read_gridbox_file(filepath):
    """
    Read gridbox parameters from a text file.
    Returns a dictionary containing the gridbox parameters.
    
    Expected format:
    PDB file: <pdb_name>
    size_x: <value>
    size_y: <value>
    size_z: <value>
    center_x: <value>
    center_y: <value>
    center_z: <value>
    """
    gridbox_params = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if ':' in line:
                    key, value = line.strip().split(':')
                    key = key.strip().lower()  # Convert to lowercase for consistent matching
                    value = value.strip()
                    if key != 'pdb file':  # Skip the PDB file line
                        # Remove any possible unit or additional text and convert to string
                        value = str(float(value.split()[0]))
                        gridbox_params[key] = value
        
        # Verify all required parameters are present
        required_params = ['size_x', 'size_y', 'size_z', 'center_x', 'center_y', 'center_z']
        for param in required_params:
            if param not in gridbox_params:
                raise ValueError(f"Missing required parameter: {param}")
        
        return gridbox_params
    except Exception as e:
        print(f"Error reading gridbox file: {str(e)}")
        return None

def get_config_files():
    """
    This function will create a AutoDock Vina config file for each PDBQT file in the directory.
    The gridbox parameters are read from a user-specified text file.
    """
    # Ask the user for the directory where the PDBQT files are located
    pwd = input("Enter the path to the 'pdbqt_files' directory: ")

    # Check if the path exists
    while not os.path.exists(pwd):
        pwd = input(
            "That path does not appear to exist.\nPlease enter the path to "
            "the directory containing the PDBQT files: "
        )
    os.chdir(pwd)

    # Ask for gridbox file path
    while True:
        gridbox_file = input("Enter the path to the gridbox parameters file: ")
        if gridbox_file.startswith('"') and gridbox_file.endswith('"'):
            gridbox_file = gridbox_file[1:-1]
        
        if not os.path.exists(gridbox_file):
            print("That path does not appear to exist.")
            continue
            
        gridbox_params = read_gridbox_file(gridbox_file)
        if gridbox_params:
            break

    # Ask the user for the remaining configuration file information
    while True:
        ligand_name = input(
            "Enter the name of the ligand PDBQT file that will be used for the docking simulation: "
        )
        # If the ligand name ends with .pdbqt, remove it
        if ligand_name.endswith(".pdbqt"):
            ligand_name = ligand_name.split(".pdbqt")[0]
        scoring = input("Enter the scoring function (ad4, vina [default] or vinardo): ")
        num_modes = input("Enter the number of modes: ")
        energy_range = input("Enter the energy range: ")
        exhaustiveness = input("Enter the exhaustiveness (default is 8; max = 32): ")
        
        # Redisplay all of this information to the user and ask them to confirm that it is correct
        print("\nIS THE FOLLOWING INFORMATION CORRECT?")
        print("Ligand name: " + ligand_name)
        print("Scoring: " + scoring)
        print("Size x: " + gridbox_params['size_x'])
        print("Size y: " + gridbox_params['size_y'])
        print("Size z: " + gridbox_params['size_z'])
        print("Center x: " + gridbox_params['center_x'])
        print("Center y: " + gridbox_params['center_y'])
        print("Center z: " + gridbox_params['center_z'])
        print("Spacing: 1")
        print("Number of modes: " + num_modes)
        print("Energy range: " + energy_range)
        print("Exhaustiveness: " + exhaustiveness)
        confirm = input("(y/n): ")
        confirm = confirm.lower()
        # Check if the user input is valid
        while confirm != "y" and confirm != "n":
            print("You did not enter a valid input.")
            confirm = input("(y/n): ")
            confirm = confirm.lower()
        if confirm == "y":
            break

    while True:
        # Create separate lists to store all of the config file information
        flex_pdbqt_files = []
        rigid_pdbqt_files = []
        pdbqt_files = [i for i in os.listdir(pwd) if i.endswith(".pdbqt") and ligand_name not in i and "residues" not in i]
        # Ask if there are any flexible residues
        flex = input("Are there any flexible residues? (y/n): ")
        flex = flex.lower()
        # Check if the user input is valid
        while flex not in ["y", "n"]:
            print("You did not enter a valid input.")
            flex = input("Are there any flexible residues? (y/n): ")
            flex = flex.lower()
        if flex == "y":
            # Prompt user for the path to the flex_residues.txt file
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
            # Read the flex_residues.txt file
            with open(residues, "r") as f:
                res_list = f.readlines()[1:]  # Remove header line

            # Loop through the files in the directory
            for file in pdbqt_files:
                # Get the base name of the file by removing "_rigid.pdbqt" or "_flex.pdbqt" or just ".pdbqt"
                base_name = file.split("_rigid.pdbqt")[0]
                base_name = base_name.split("_flex.pdbqt")[0]
                base_name = base_name.split(".pdbqt")[0]

                for line in res_list:
                    if line.split(":")[0] == base_name:
                        res = line.split("[")[1].split("]")[0]
                        if res == "":
                            rigid_pdbqt_files.append(base_name)
                        else:
                            flex_pdbqt_files.append(base_name)
        else:
            # If there are no flexible residues, add all to the rigid_pdbqt_files list
            for file in pdbqt_files:
                base_name = file.split(".pdbqt")[0]
                rigid_pdbqt_files.append(base_name)
                
        # Remove duplicates in the lists
        flex_pdbqt_files = list(set(flex_pdbqt_files))
        rigid_pdbqt_files = list(set(rigid_pdbqt_files))

        # Ask the user if the entries are correct
        print("\nIS THE FOLLOWING INFORMATION CORRECT?")
        print("FLEXIBLE PDBQTs: ")
        for i in flex_pdbqt_files:
            print(i)
        print("RIGID PDBQTs: ")
        for i in rigid_pdbqt_files:
            print(i)
        confirm = input("(y/n): ")
        confirm = confirm.lower()
        # Check if the user input is valid
        while confirm != "y" and confirm != "n":
            print("You did not enter a valid input.")
            confirm = input("(y/n): ")
            confirm = confirm.lower()
        if confirm == "y":
            break
    
    # Modified file writing section for flexible PDBQT files
    for file in flex_pdbqt_files:
        with open(file.split("_", 1)[0] + "_conf" + ".txt", "w") as f:
            f.write("flex = {}_flex.pdbqt\n".format(file))
            f.write("receptor = {}_rigid.pdbqt\n".format(file))
            f.write("ligand = {}.pdbqt\n".format(ligand_name))
            f.write("scoring = {}\n\n".format(scoring))
            f.write("size_x = {}\n".format(gridbox_params['size_x']))
            f.write("size_y = {}\n".format(gridbox_params['size_y']))
            f.write("size_z = {}\n\n".format(gridbox_params['size_z']))
            f.write("center_x = {}\n".format(gridbox_params['center_x']))
            f.write("center_y = {}\n".format(gridbox_params['center_y']))
            f.write("center_z = {}\n\n".format(gridbox_params['center_z']))
            f.write("spacing = 1\n\n")
            f.write("num_modes = {}\n".format(num_modes))
            f.write("energy_range = {}\n".format(energy_range))
            f.write("exhaustiveness = {}\n\n".format(exhaustiveness))
            f.write("out = " + file.split("_", 1)[0] + "_bound_" + ligand_name + ".pdbqt")

    # Modified file writing section for rigid PDBQT files
    for file in rigid_pdbqt_files:
        with open(file.split("_", 1)[0] + "_conf" + ".txt", "w") as f:
            f.write("receptor = {}.pdbqt\n".format(file))
            f.write("ligand = {}.pdbqt\n".format(ligand_name))
            f.write("scoring = {}\n\n".format(scoring))
            f.write("size_x = {}\n".format(gridbox_params['size_x']))
            f.write("size_y = {}\n".format(gridbox_params['size_y']))
            f.write("size_z = {}\n\n".format(gridbox_params['size_z']))
            f.write("center_x = {}\n".format(gridbox_params['center_x']))
            f.write("center_y = {}\n".format(gridbox_params['center_y']))
            f.write("center_z = {}\n\n".format(gridbox_params['center_z']))
            f.write("spacing = 1\n\n")
            f.write("num_modes = {}\n".format(num_modes))
            f.write("energy_range = {}\n".format(energy_range))
            f.write("exhaustiveness = {}\n\n".format(exhaustiveness))
            f.write("out = " + file.split("_", 1)[0] + "_bound_" + ligand_name + ".pdbqt")

    print("Configuration files have been saved to the 'pdbqt_files' directory.")

if __name__ == "__main__":
    get_config_files()