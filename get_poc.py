"""
 This script allows you to get the protein binding pocket residues for a 
 given ligand in a PDB file. The methods are based on the PLIP program                                              
 (https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index). The script will 
 calculate the distance between the ligand and each residue in the PDB file and 
 return the residues within a specified distance of the ligand. The default distance 
 is 7.5 Angstroms. The PDB file should NOT be protonated!                                                           
 This script will return a binding pocket PDB file and .poc file.                                 
"""

import os
import shutil

import numpy as np


def euclidean3d(v1, v2):
    """
    Faster implementation of euclidean distance for the 3D case.
    """
    if not len(v1) == 3 and len(v2) == 3:
        return None
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def centroid(coo):
    """
    Calculates the centroid from a 3D point cloud and returns the coordinates
    :param coo: Array of coordinate arrays
    :returns : centroid coordinates as list
    """
    return list(
        map(
            np.mean,
            (([c[0] for c in coo]), ([c[1] for c in coo]), ([c[2] for c in coo])),
        )
    )


def min_dist(pdb_file, ligand):
    """
    Determines the residues within 7.5 Angstroms of the ligand in the PDB file.
    """
    # Create a list to hold the coordinates of the ligand atoms
    ligand_coords = []
    # Read PDB file
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        # Add the coordinates of the ligand atoms to the ligand_coords list
        if line.startswith("HETATM"):
            if line[17:20].strip() == ligand:
                ligand_coords.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
    # Create a dictionary to store the residue coordinates
    res_coords = {}
    # Iterate through the lines of the PDB file on a per residue basis to calculate 
    # the residue centroids
    for line in lines:
        if line.startswith("ATOM"):
            # Extract the residue number (resnr) from the line
            resnr = line[22:26]
            # Extract the coordinates from the line
            coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            # If the resnr is not in the res_coords dictionary, add it and set the 
            # value to a list containing the coordinates
            if not resnr in res_coords:
                res_coords[resnr] = [coords]
            # If the resnr is in the res_coords dict, append the coordinates to the list
            else:
                res_coords[resnr].append(coords)
    # Calculate the minimum distance between the each of the ligand and the residue atoms 
    # on a per residue basis
    min_dist = {}
    for resnr, coords in res_coords.items():
        min_dist[resnr] = 100
        for ligand_coord in ligand_coords:
            for coord in coords:
                dist = euclidean3d(ligand_coord, coord)
                if dist < min_dist[resnr]:
                    min_dist[resnr] = dist
    # Save only the residue numbers that have a min distance less than 7.5 Angstroms
    min_dist = {resnr: dist for resnr, dist in min_dist.items() if dist < 7.5}
    return min_dist


def get_bindingsite(pdb_file, ligand):
    """
    Determines which residues in the PDB model are within 7.5 Angstroms of the ligand
    according to an all-atom distance calculation. Additionally, the centroid of the 
    ligand and the centroid of each residue are calculated. The distance between the 
    ligand centroid and the residue centroid is calculated and the residues with a 
    distance less than 7.5 Angstroms plus the ligand's radius are returned.
    """
    # Create a list to hold the coordinates of the ligand atoms
    ligand_coords = []
    # Read PDB file
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        # Add the coordinates of the ligand atoms to the ligand_coords list
        if line.startswith("HETATM"):
            if line[17:20].strip() == ligand:
                ligand_coords.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
    # Check if ligand coordinates were found
    if not ligand_coords:
        print(f"Ligand {ligand} not found in PDB file {pdb_file}.")
        return None, None
    # Calculate the centroid of the ligand using the ligand coordinates 
    # and the centroid function
    ligand_centroid = centroid(ligand_coords)
    print(f"Ligand centroid: {ligand_centroid}")
    # Create a dictionary to store the residue coordinates
    res_coords = {}
    # Iterate through the lines of the PDB file on a per residue basis to 
    # calculate the residue centroids
    for line in lines:
        if line.startswith("ATOM"):
            # Extract the residue number (resnr) from the line
            resnr = line[22:26]
            # Extract the coordinates from the line
            coords = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            # If the resnr is not in the res_centroids dictionary, add it and 
            # set the value to a list containing the coordinates
            if not resnr in res_coords:
                res_coords[resnr] = [coords]
            # If the resnr is in the res_centroids dictionary, append the 
            # coordinates to the list
            else:
                res_coords[resnr].append(coords)
    # Create a dictionary to store the residue centroids
    res_centroids = {}
    # Iterate through the res_coords dictionary and calculate the 
    # centroid of each residue
    for resnr, coords in res_coords.items():
        res_centroids[resnr] = centroid(coords)
    # Define the binding site distance (BS_DIST) in Angstroms
    BS_DIST = 7.5
    # Calculate the distance between the ligand centroid and the ligand 
    # atoms and store the maximum distance in max_dist
    max_dist = 0
    for coords in ligand_coords:
        dist = euclidean3d(ligand_centroid, coords)
        if dist > max_dist:
            max_dist = dist
    # Define the cutoff distance as the BS_DIST + max_dist
    cutoff = BS_DIST + max_dist
    # Print the cutoff distance
    print(f"Cutoff distance: {cutoff}")
    # Create a list to store the binding site residue numbers
    bindingsite_resnr = []
    # Iterate through the res_centroids and add the resnr to the 
    # bindingsite_resnr list if the distance between the ligand centroid 
    # and the residue centroid is less than or equal to the cutoff
    for resnr, coords in res_centroids.items():
        if euclidean3d(ligand_centroid, coords) < cutoff:
            bindingsite_resnr.append(resnr)
    # Remove duplicate resnr from the bindingsite_resnr list
    bindingsite_resnr = list(set(bindingsite_resnr))
    # Calculate the minimum distance between the ligand and each 
    # residue in the binding site
    min_dist_res = min_dist(pdb_file, ligand)
    # Remove resnr that are not matching with the min_dist_res list 
    # and the current bindingsite_resnr list
    bindingsite_resnr = [i for i in bindingsite_resnr if i in min_dist_res]
    # Sort the bindingsite_resnr list
    bindingsite_resnr.sort()
    # Count the number of residues in the binding site
    print(f"Number of residues in the binding site: {len(bindingsite_resnr)}")
    # Return the bindingsite_resnr list
    return bindingsite_resnr


def get_pdb_bindingsite(pdb_file, bindingsite_res):
    """
    Extracts the binding site residues from the PDB file and returns
    a new PDB file containing only the binding site residues.
    """
    # Convert resnr to int
    bindingsite_res = [int(i) for i in bindingsite_res]
    # Read PDB file
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    # Create a list to store the new POC lines
    poc = []
    # Iterate through the lines of the PDB file
    for line in lines:
        # Check if the line is an ATOM or HETATM record
        if line.startswith("ATOM"):
            # Extract the residue number (resnr) from the line
            resnr = line[22:26]
            # convert resnr to an integer
            resnr = int(resnr)
            # If this resnr is equivalent to any resnr in the 
            # binding_site_resnr list, add the line to the new PDB list
            if resnr in bindingsite_res:
                poc.append(line)
        # Check if the line starts with 'TER' and add it to the new PDB list
        if line.startswith("TER"):
            # Remove all characters after the first 3
            line = line[:3]
            poc.append(line)
        # Skip all other lines
        else:
            continue
    # Get the working directory of the original PDB file
    working_dir = os.path.dirname(pdb_file)
    # Get the filename of the original PDB file without the extension
    filename = os.path.basename(pdb_file).split(".")[0]
    # Add a title line to the new PDB list
    poc.insert(0, "POC  " + filename + "_poc\n")
    # Create a new filename for the new PDB file
    poc_filename = filename.split(".")[0] + "_bindingsite.poc"
    # Create a new PDB file in the working directory of the original PDB file
    with open(os.path.join(working_dir, poc_filename), "w") as f:
        for line in poc:
            f.write(line)
    # Create a new filename for the new PDB file
    pdb_filename = filename.split(".")[0] + "_bindingsite.pdb"
    # Create a new PDB file in the working directory of the original PDB file
    with open(os.path.join(working_dir, pdb_filename), "w") as f:
        for line in poc:
            f.write(line)
    # Print the path to the new PDB file
    print(f"New PDB file: {os.path.join(working_dir, pdb_filename)}")

def get_gridbox_coordinates(pdb_file, ligand):
    """
    Calculates the gridbox size and coordinates for a given ligand in a PDB file.
    It does this by calculating the maximum distance between the ligand centroid
    and the ligand atoms. The gridbox size is then calculated as 4 times the maximum
    distance in each dimension (x, y, z). The gridbox coordinates are then calculated
    as the ligand centroid coordinates.
    """
    # Create a list to hold the coordinates of the ligand atoms
    ligand_coords = []
    # Read PDB file
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        # Add the coordinates of the ligand atoms to the ligand_coords list
        if line.startswith("HETATM"):
            if line[17:20].strip() == ligand:
                ligand_coords.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
    # Check if ligand coordinates were found
    if not ligand_coords:
        print(f"The ligand ({ligand}) was not found in the given PDB file {pdb_file}.")
        return None, None
    # Calculate the centroid of the ligand using the ligand coordinates 
    # and the centroid function
    ligand_centroid = centroid(ligand_coords)
    # Round the ligand centroid coordinates to 3 decimal places
    ligand_centroid = [round(i, 3) for i in ligand_centroid]
    # Calculate the maximum distance between the ligand centroid and the ligand atoms
    max_dist = 0
    for coords in ligand_coords:
        dist = euclidean3d(ligand_centroid, coords)
        if dist > max_dist:
            max_dist = dist
    # Using the max_dist, calculate the size of the gridbox
    size = round(max_dist * 4, 3)
    # Write the ligand centroid coordinates to a file
    with open(pdb_file.split(".")[0] + "_gridbox_coords.txt", "w") as f:
        f.write(
            "PDB file: "
            + pdb_file.split(".")[0]
            + "\n\nsize_x: "
            + str(size)
            + "\nsize_y: "
            + str(size)
            + "\nsize_z: "
            + str(size)
            + "\n\ncenter_x: "
            + str(ligand_centroid[0])
            + "\ncenter_y: "
            + str(ligand_centroid[1])
            + "\ncenter_z: "
            + str(ligand_centroid[2])
        )


def check_multi_lig(pdb_file, ligand):
    """
    Checks if there are multiple ligands in a PDB file.
    """
    # Read PDB file
    with open(pdb_file, "r") as f:
        lines = f.readlines()
    # Create a list to store the ligand atom names
    ligand_atoms = []
    # Iterate through the lines of the PDB file
    for line in lines:
        # Add the ligand atom names to the ligand_atoms list
        if line.startswith("HETATM"):
            if line[17:20].strip() == ligand:
                atom = line[12:16].strip()
                chain = line[21].strip()
                if atom in ligand_atoms:
                    return True
                else:
                    ligand_atoms.append(atom)
    return False


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


def get_poc():
    """
    Main function
    """
    # Ask the user if they would like to analyze a single PDB file 
    # or a directory of PDB files
    pdb_type = input(
        "Analyze a single PDB file or a directory of PDB files? (single/directory): "
    )
    if check_exit(pdb_type):
        return
    # Check if the user made a typo and if so, ask them to enter the input again
    while pdb_type not in ["single", "directory"]:
        print('Invalid input. Please enter either "single" or "directory".')
        pdb_type = input(
            "Analyze a single PDB file or a directory of PDB files? (single/directory): "
        )
    # Ask the user to enter the file directory containing the PDB files
    pdb_dir = input("Enter the file directory containing the PDB file(s): ")
    if check_exit(pdb_dir):
        return
    # Change the working directory to the directory containing the PDB files
    os.chdir(pdb_dir)
    # If the user would like to analyze a single PDB file
    if pdb_type == "single":
        # Ask the user to enter the PDB file
        pdb_file = input("Enter the PDB file: ")
        # Check if the user included the file extension
        if not pdb_file.endswith(".pdb"):
            pdb_file += ".pdb"
        # Check if the PDB file exists and if not, ask the user to enter the PDB file again
        while not os.path.exists(pdb_file):
            print(
                f"{pdb_file} does not appear to exist. Please enter a valid PDB file."
                " The entry is case sensitive."
            )
            pdb_file = input("Enter the PDB file: ")
            if not pdb_file.endswith(".pdb"):
                pdb_file += ".pdb"
            if check_exit(pdb_file):
                return
                # Ask the user to enter the ligand name
        ligand = input("Enter the ligand ID as it is found within the PDB: ")
        if check_exit(ligand):
            return
        # Read the PDB file
        with open(pdb_file, "r") as f:
            lines = f.readlines()
        # Check if the ligand is present in the PDB file, 
        # if not ask the user to enter the ligand name again
        while ligand not in [
            line[17:20].strip() for line in lines if line.startswith("HETATM")
        ]:
            print(
                f"{os.path.basename(pdb_file)} does not appear to contain "
                "the ligand {ligand}. Please try entering it again."
            )
            # Try to get the ligand name again
            ligand = input("Enter the ligand ID as it is found within the PDB: ")
            if check_exit(ligand):
                return
                # Check if there are multiple of the same ligand in the PDB file
        if check_multi_lig(pdb_file, ligand):
            print(
                f"ERROR: It appears that there are multiple {ligand} molecules "
                "in {pdb_file}.\nPlease remove all but one of the {ligand} molecules "
                "from the PDB file before using this function."
            )
            return
        # Check if the PDB file is protonated
        with open(pdb_file, "r") as f:
            lines = f.readlines()
        if any(
            [
                line[12:16].strip().startswith("H")
                for line in lines
                if line.startswith("ATOM")
            ]
        ):
            print(
                f"It appears that {pdb_file} is protonated.\nThe best results will "
                "be obtained by using an unprotonated PDB file."
            )
            # Ask the user if they would like to continue or use a different PDB file
            while True:
                continue_pdb = input("Continue with the current PDB file? (y/n): ")
                if continue_pdb == "y":
                    break
                elif continue_pdb == "n":
                    pdb_file = input("Enter the PDB file: ")
                    if check_exit(pdb_file):
                        return
                    if not pdb_file.endswith(".pdb"):
                        pdb_file += ".pdb"
                        while not os.path.exists(pdb_file):
                            print(
                                f"{pdb_file} does not exist. Please enter a "
                                "valid PDB file."
                            )
                            pdb_file = input("Enter the PDB file: ")
                            if check_exit(pdb_file):
                                return
                elif check_exit(continue_pdb):
                    return
                else:
                    print('Invalid input. Please enter either "y" or "n".')
        # Get the binding site residues
        bindingsite_res = get_bindingsite(pdb_file, ligand)
        # Generate a new PDB file with the binding site residues
        get_pdb_bindingsite(pdb_file, bindingsite_res)
        # Retain the 'HETATM' records in the new PDB file?
        retain_hetatm = input(
            "Retain the ligand records in the new PDB binding pocket file? (y/n): "
        )
        # Check if the user made a typo and if so, ask them to enter the input again
        while retain_hetatm not in ["y", "n"]:
            print('Invalid input. Please enter either "y" or "n".')
            retain_hetatm = input(
                "Retain the ligand records in the new PDB file? (y/n): "
            )
            if check_exit(retain_hetatm):
                return
        # If the user would like to retain the ligand records, transfer the ligand 
        # records from the original PDB file to the new PDB file
        if retain_hetatm == "y":
            # Read the original PDB file
            with open(pdb_file, "r") as f:
                lines = f.readlines()
                # Get the ligand records
                hetatm_records = [
                    line
                    for line in lines
                    if line.startswith("HETATM") and line[17:20].strip() == ligand
                ]
                # Add the ligand records to the new PDB file
                with open(
                    os.path.join(
                        os.path.dirname(pdb_file),
                        os.path.basename(pdb_file).split(".")[0] + "_bindingsite.pdb",
                    ),
                    "a",
                ) as f:
                    f.write("\n")
                    for line in hetatm_records:
                        f.write(line)
            print(
                "The ligand records have been retained in the new "
                "binding site PDB file."
            )
        else:
            print(
                "The ligand records have been removed from the new "
                "binding site PDB file."
            )
        # Read the new PDB file
        with open(
            os.path.join(
                os.path.dirname(pdb_file),
                os.path.basename(pdb_file).split(".")[0] + "_bindingsite.pdb",
            ),
            "r",
        ) as f:
            lines = f.readlines()
            # Move the 'TER' record to the end of the file
            for line in lines:
                if line.startswith("TER"):
                    lines.remove(line)
                    lines.append(line)
        # Write the new PDB file
        with open(
            os.path.join(
                os.path.dirname(pdb_file),
                os.path.basename(pdb_file).split(".")[0] + "_bindingsite.pdb",
            ),
            "w",
        ) as f:
            for line in lines:
                f.write(line)
        # Ask the user if they would like to output the gridbox coordinates
        output_gridbox = input("Output the gridbox coordinates? (y/n): ")
        if check_exit(output_gridbox):
            return
        # Check if the user made a typo and if so, ask them to enter the input again
        while output_gridbox not in ["y", "n"]:
            print('Invalid input. Please enter either "y" or "n".')
            output_gridbox = input("Output the gridbox coordinates? (y/n): ")
            if check_exit(output_gridbox):
                return
        # If the user would like to output the gridbox coordinates, 
        # run get_gridbox_coordinates()
        if output_gridbox == "y":
            get_gridbox_coordinates(pdb_file, ligand)
            # If there is a "details" folder, move the "bindingsite" PDBs and the gridbox 
            # coordinates to the "details" folder
            if os.path.exists(os.path.join(os.path.dirname(pdb_file), "details")):
                shutil.move(
                    os.path.join(
                        os.path.dirname(pdb_file),
                        os.path.basename(pdb_file).split(".")[0] + "_bindingsite.pdb",
                    ),
                    os.path.join(os.path.dirname(pdb_file), "details"),
                )
                shutil.move(
                    os.path.join(
                        os.path.dirname(pdb_file),
                        os.path.basename(pdb_file).split(".")[0] + "_bindingsite.poc",
                    ),
                    os.path.join(os.path.dirname(pdb_file), "details"),
                )
                shutil.move(
                    os.path.join(
                        os.path.dirname(pdb_file),
                        os.path.basename(pdb_file).split(".")[0]
                        + "_gridbox_coords.txt",
                    ),
                    os.path.join(os.path.dirname(pdb_file), "details"),
                )
                print(
                    'The gridbox coordinates and bindingsite files have been '
                    'outputted to the "details" folder.'
                )
            else:
                print(
                    "The gridbox coordinates and bindingsite files have been "
                    "outputted to the current directory."
                )
    # If the user would like to analyze a directory of PDB files
    elif pdb_type == "directory":
        # Ask the user to enter the ligand name
        ligand = input("Enter the ligand name: ")
        # Create a new folder in the pdb_dir to store the binding site PDB files
        if not os.path.exists(os.path.join(pdb_dir, "binding_sites")):
            os.mkdir(os.path.join(pdb_dir, "binding_sites"))
        # Get the list of PDB files in the directory
        pdb_files = [
            os.path.join(pdb_dir, i) for i in os.listdir(pdb_dir) if i.endswith(".pdb")
        ]
        # Iterate through the PDB files and generate a new PDB file 
        # with the binding site residues
        for pdb_file in pdb_files:
            # Check if the file name contains "bindingsite", if so, skip it
            if "bindingsite" in os.path.basename(pdb_file):
                continue
            # If the PDB file does not contain the ligand, skip to the next PDB file
            with open(pdb_file, "r") as f:
                lines = f.readlines()
            # Inform the user which PDB file is being processed, 
            # and extract the PDB file name.
            print(f"Processing {os.path.basename(pdb_file)}")
            # Check if the ligand is present in the PDB file
            if ligand not in [
                line[17:20].strip() for line in lines if line.startswith("HETATM")
            ]:
                print(
                    f"{os.path.basename(pdb_file)} does not contain the "
                    "ligand {ligand}. Skipping to the next PDB file."
                )
                continue
            # Check if the PDB file is protonated
            if any(
                [
                    line[12:16].strip().startswith("H")
                    for line in lines
                    if line.startswith("ATOM")
                ]
            ):
                print(
                    f"It appears that {os.path.basename(pdb_file)} is "
                    "protonated.\nSkipping to the next PDB file."
                )
                # Skip the file if it is protonated
                continue
            # Get the binding site residue IDs
            bindingsite_res = get_bindingsite(pdb_file, ligand)
            # Generate a new PDB file with the binding site residues
            get_pdb_bindingsite(pdb_file, bindingsite_res)
            # Move the new PDB and POC files to the binding_sites folder
            shutil.move(
                os.path.join(
                    os.path.dirname(pdb_file),
                    os.path.basename(pdb_file).split(".")[0] + "_bindingsite.pdb",
                ),
                os.path.join(pdb_dir, "binding_sites"),
            )
            shutil.move(
                os.path.join(
                    os.path.dirname(pdb_file),
                    os.path.basename(pdb_file).split(".")[0] + "_bindingsite.poc",
                ),
                os.path.join(pdb_dir, "binding_sites"),
            )
        # Inform the user where the new files were stored
        print(
            f'\nThe binding site PDB and POC files have been '
            'stored in the "binding_sites" folder.'
        )
        # Combine the .poc files into one file?
        combine = input(
            "\nWould you like to combine the query .poc files into one file? (y/n): "
        )
        while combine not in ["y", "n"]:
            print('Invalid input. Please enter either "y" or "n".')
            combine = input("\nWould you like to combine the query .poc files into one file? (y/n): ")
            if check_exit(combine):
                return
        if combine == "y":
            os.chdir(os.path.join(pdb_dir, "binding_sites"))
            # Create a new empty .poc file
            ref_poc = input(
                "\nEnter the name of the reference .poc file so that it "
                "doesn't get added to the combined file: "
            )
            if not ref_poc.endswith(".pdb"):
                ref_poc += ".pdb"
            # Check if the reference .poc file exists
            while not os.path.exists(os.path.join(pdb_dir, ref_poc)):
                print(
                    f'The file "{ref_poc}" does not exist. Please enter the '
                    "name of the reference .poc file so that it doesn't get "
                    "added to the combined file: "
                )
                ref_poc = input()
                if not ref_poc.endswith(".pdb"):
                    ref_poc += ".pdb"
                if check_exit(ref_poc):
                    return
            all_poc = []
            # Read all of the .poc files and add them to the all_poc list
            for file in os.listdir(os.path.join(pdb_dir, "binding_sites")):
                if file.endswith(".poc") and file.split(".")[0] != ref_poc.split(".")[0]:
                    with open(file, "r") as f:
                        for line in f:
                            if line.startswith("TER"):
                                line += "\n"
                            all_poc.append(line)
            with open("all_pocs.poc", "w") as o:
                for line in all_poc:
                    o.write(line)
            print(
                '\nAll of the query .poc files have been combined and saved '
                'as "all_pocs.poc" in the "binding_sites" folder.'
            )
        else:
            print(
                "\nAll .poc files have been stored in the binding_sites folder."
            )


if __name__ == "__main__":
    get_poc()
