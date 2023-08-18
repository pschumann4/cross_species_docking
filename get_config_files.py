"""
 This script will create a AutoDock Vina config file for each PDBQT file in the directory.        
 The user will be prompted to enter the following information:                                    
                                                                                                  
 1. The name of the ligand PDBQT file that will be used for the docking simulation                
 2. The scoring function (ad4, vina [default] or vinardo)                                         
 3. The size of the gridbox in the x, y and z dimensions                                          
 4. The center coordinate of the gridbox in the x, y and z dimensions                             
 5. The spacing between grid points (default is 0.375)                                            
 6. The number of modes                                                                           
 7. The energy range                                                                              
 8. The exhaustiveness (default is 8; max = 32)                                                   
 9. Which PDBQT files should be treated as flexible                                               
                                                                                                  
 The user will be asked to confirm that all of the information is correct before the config files 
 are created.                                                                                     
"""

import os


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


def get_config_files():
    """
    This function will create a AutoDock Vina config file for each PDBQT file in the directory.
    The user will be prompted to enter the following information:
    1. The name of the ligand PDBQT file that will be used for the docking simulation
    2. The scoring function (ad4, vina [default] or vinardo)
    3. The size of the gridbox in the x, y and z dimensions
    4. The center coordinate of the gridbox in the x, y and z dimensions
    5. The spacing between grid points (default is 0.375)
    6. The number of modes
    7. The energy range
    8. The exhaustiveness (default is 8; max = 32)
    9. Which PDBQT files should be treated as flexible
    The user will be asked to confirm that all of the information is correct before the config files
    are created.
    """
    # Ask the user for the directory where the PDBQT files are located
    pwd = input("Enter the directory where the PDBQT structure files are located: ")
    if check_exit(pwd):
        return
    # Check that the directory exists, if not, ask the user to enter a new directory
    while not os.path.isdir(pwd):
        print("The directory you entered does not exist.")
        pwd = input("Enter the directory where the PDBQT structure files are located: ")
        if check_exit(pwd):
            return
    # Change the working directory to the directory where the PDBQT files are located
    os.chdir(pwd)
    # Ask the user for all of the configuration file information.
    while True:
        ligand_name = input(
            "Enter the name of the ligand PDBQT file that will be used for the docking simulation: "
        )
        # If the ligand name ends with .pdbqt, remove it
        if ligand_name.endswith(".pdbqt"):
            ligand_name = ligand_name.split(".pdbqt")[0]
        scoring = input("Enter the scoring function (ad4, vina [default] or vinardo): ")
        size_x = input("Enter the size of the gridbox in the x dimension: ")
        size_y = input("Enter the size of the gridbox in the y dimension: ")
        size_z = input("Enter the size of the gridbox in the z dimension: ")
        center_x = input(
            "Enter the center coordinate of the gridbox in the x dimension: "
        )
        center_y = input(
            "Enter the center coordinate of the gridbox in the y dimension: "
        )
        center_z = input(
            "Enter the center coordinate of the gridbox in the z dimension: "
        )
        spacing = input("Enter the spacing between grid points (default is 0.375): ")
        num_modes = input("Enter the number of modes: ")
        energy_range = input("Enter the energy range: ")
        exhaustiveness = input("Enter the exhaustiveness (default is 8; max = 32): ")
        # Redisplay all of this information to the user and ask them to confirm that it is correct
        print("\nIS THE FOLLOWING INFORMATION CORRECT?")
        print("Ligand name: " + ligand_name)
        print("Scoring: " + scoring)
        print("Size x: " + size_x)
        print("Size y: " + size_y)
        print("Size z: " + size_z)
        print("Center x: " + center_x)
        print("Center y: " + center_y)
        print("Center z: " + center_z)
        print("Spacing: " + spacing)
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
    # Get flexibility info
    while True:
        # Create separate lists to store all of the config file information
        flex_pdbqt_files = []
        rigid_pdbqt_files = []
        pdbqt_files = [i for i in os.listdir(pwd) if i.endswith(".pdbqt") and ligand_name not in i and "residues" not in i]
        # Loop through the files in the directory
        for file in pdbqt_files:
            # Get the base name of the file
            base_name = file.split(".pdbqt")[0]
            # Ask the user if the file is flexible or rigid
            is_flex = input(
                "Should " + base_name + " be treated as flexible? (y/n): "
            )
            # Convert the user input to lower case
            is_flex = is_flex.lower()
            # Check if the user input is valid
            while is_flex != "y" and is_flex != "n":
                print("You did not enter a valid input.")
                is_flex = input(
                    "Should " + base_name + " be treated as flexible? (y/n): "
                )
                is_flex = is_flex.lower()
            if is_flex == "y":
                # Add the base name to the flex_pdbqt_files list
                flex_pdbqt_files.append(base_name)
            elif is_flex == "n":
                # Add the base name to the rigid_pdbqt_files list
                rigid_pdbqt_files.append(base_name)
        # Ask the user if the entries are correct
        print("\nIS THE FOLLOWING INFORMATION CORRECT?")
        print("Flexible PDBQT files: " + str(flex_pdbqt_files))
        print("Rigid PDBQT files: " + str(rigid_pdbqt_files))
        confirm = input("(y/n): ")
        confirm = confirm.lower()
        # Check if the user input is valid
        while confirm != "y" and confirm != "n":
            print("You did not enter a valid input.")
            confirm = input("(y/n): ")
            confirm = confirm.lower()
        if confirm == "y":
            break
    # Loop through the flexible PDBQT files and create a config file for each one
    for file in flex_pdbqt_files:
        with open(file.split("_", 1)[0] + "_conf" + ".txt", "w") as f:
            f.write("flex = {}_flex_residues.pdbqt\n".format(file))
            f.write("receptor = {}_rigid_residues.pdbqt\n".format(file))
            f.write("ligand = {}.pdbqt\n".format(ligand_name))
            f.write("scoring = {}\n\n".format(scoring))
            f.write("size_x = {}\n".format(size_x))
            f.write("size_y = {}\n".format(size_y))
            f.write("size_z = {}\n\n".format(size_z))
            f.write("center_x = {}\n".format(center_x))
            f.write("center_y = {}\n".format(center_y))
            f.write("center_z = {}\n\n".format(center_z))
            f.write("spacing = {}\n\n".format(spacing))
            f.write("num_modes = {}\n".format(num_modes))
            f.write("energy_range = {}\n".format(energy_range))
            f.write("exhaustiveness = {}\n\n".format(exhaustiveness))
            f.write("out = " + file.split("_", 1)[0] + "_bound_" + ligand_name + ".pdbqt")
    # Loop through the rigid PDBQT files and create a config file for each one
    for file in rigid_pdbqt_files:
        with open(file.split("_", 1)[0] + "_conf" + ".txt", "w") as f:
            f.write("receptor = {}.pdbqt\n".format(file))
            f.write("ligand = {}.pdbqt\n".format(ligand_name))
            f.write("scoring = {}\n\n".format(scoring))
            f.write("size_x = {}\n".format(size_x))
            f.write("size_y = {}\n".format(size_y))
            f.write("size_z = {}\n\n".format(size_z))
            f.write("center_x = {}\n".format(center_x))
            f.write("center_y = {}\n".format(center_y))
            f.write("center_z = {}\n\n".format(center_z))
            f.write("spacing = {}\n\n".format(spacing))
            f.write("num_modes = {}\n".format(num_modes))
            f.write("energy_range = {}\n".format(energy_range))
            f.write("exhaustiveness = {}\n\n".format(exhaustiveness))
            f.write("out = " + file.split("_", 1)[0] + "_bound_" + ligand_name + ".pdbqt")
    # Check that there is a config file for each PDBQT file
    for file in os.listdir(pwd):
        if (
            file.endswith(".pdbqt")
            and ligand_name not in file
            and "residues" not in file
        ):
            base_name = file.split(".pdbqt")[0].split("_", 1)[0]
            if base_name + "_conf" + ".txt" not in os.listdir(pwd):
                print("There is no config file for " + base_name + ".pdbqt")
            else:
                print("Configuration file was created for " + base_name)


if __name__ == "__main__":
    get_config_files()
