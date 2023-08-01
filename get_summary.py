####################################################################################################
# This script will create a summary file containing the binding affinity, PPS-Score, ligand RMSD,  #
# and PLIF Tanimoto values for each binding model. The user will be prompted to enter the          #
# following information:                                                                           #
#                                                                                                  #
# 1. The directory containing the AutoDock Vina log files                                          #
# 2. The name of the ligand, which should match the suffix of the log files                        #
# 3. The file path for the PPS-Score file                                                          #
# 4. The file path for the ligand RMSD file                                                        #
# 5. The file path for the PLIF Tanimoto matrix file                                               #
# 6. The name of the reference model (as listed in the matrix)                                     #
# 7. The output directory                                                                          #
#                                                                                                  #
# The summary file will be written to the specified output directory.                              #
####################################################################################################

############################
# Import required packages #
############################

import os
import pandas as pd

#################
# Load function #
#################

def get_summary():
    # Initialize empty dataframe to store the summary data
    summary_df = pd.DataFrame(columns=["binding_model", "species", "pose", "binding_affinity", "ppsscore", "lig_rmsd", "plif_tanimoto"])
    # Prompt user for directory containing the AutoDock Vina log files
    vina_logs = input("Enter the file directory contating the AutoDock Vina log files: ")
    # Prompt user for ligand name
    ligand = input("Enter the name of the ligand, which should match the suffix of the log files: ")
    # Initialize lists to hold binding affinities, poses, and species names
    binding_affinities = []
    poses = []
    species = []
    file_count = 0

    # Read the AutoDock Vina log files
    for file in os.listdir(vina_logs):
        # If the file contains the words "bound" and "ligand" and endswith ".pdbqt"
        if file.endswith(".pdbqt") and "bound" in file and "ligand" in file:
            file_count += 1
            with open(os.path.join(vina_logs, file), "r") as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("REMARK VINA RESULT:"):
                        # Get the binding affinity
                        binding_affinity = line.split()[3]
                        # Append the binding affinity and pose to the lists
                        binding_affinities.append(binding_affinity)
                        pose = file.split("_")[-1].split(".")[0]
                        poses.append(pose)
                        # Get the species name by splitting the file name at the underscore
                        species_name = file.split("_")[0]
                        # Append the species name to the list
                        species.append(species_name)

    # Add the binding affinities, poses, and species names to the dataframe
    summary_df["binding_affinity"] = binding_affinities
    summary_df["pose"] = poses
    summary_df["species"] = species

    # Get the number of rows in the dataframe
    n_rows = len(summary_df.index)
    binding_models = list(range(1, n_rows + 1))
    # Add the binding models to the dataframe
    summary_df["binding_model"] = binding_models
    # Print the df
    print(summary_df)

    # Prompt user for the directory containing the PPS-Score files
    ppsscore_file = input("Enter the file path for the PPS-Score file: ")
    ppsscore_file = ppsscore_file.replace('"', "")
    # Create a list to hold the PPS-Score values
    ppsscores = []

    # Read the PPS-Score file
    with open(ppsscore_file, "r") as f:
        lines = f.readlines()
        line_counter = 2
        for line in lines[2:file_count + 2]:
            pps_line = lines[line_counter]
            pps_line = pps_line.split()
            if len(pps_line) != 0:
                ppsscores.append(pps_line[2])
            if line_counter < file_count + 1:
                line_counter += 1

    # Add the PPS-Score values to the dataframe
    summary_df["ppsscore"] = ppsscores
    # Print the df
    print(summary_df)

    # Prompt user for the ligand RMSD file
    lig_rmsd_file = input("Enter the file path for the ligand RMSD file: ")
    # Remove quotes from the file path
    lig_rmsd_file = lig_rmsd_file.replace('"', "")
    # Read the ligand RMSD file
    with open(lig_rmsd_file, "r") as f:
        lines = f.readlines()
        lines = [line.split() for line in lines]
        lig_rmsds = [line[1] for line in lines[2:]]
    # Add the ligand RMSD values to the dataframe
    summary_df["lig_rmsd"] = lig_rmsds
    # Print the df
    print(summary_df)

    # Prompt user for the PLIF Tanimoto file
    plif_tanimoto_file = input("Enter the file path for the PLIF Tanimoto matrix file: ")
    plif_tanimoto_file = plif_tanimoto_file.replace('"', "")
    plif_tanimoto_df = pd.read_csv(plif_tanimoto_file)
    # Prompt user for the name of the reference model
    ref_model = input("Enter the name of the reference model (as listed in the matrix): ")
    # Remove the column corresponding to the reference model
    plif_tanimoto_df = plif_tanimoto_df.drop(ref_model, axis=1)
    plif_tanimoto_df = plif_tanimoto_df.set_index(plif_tanimoto_df.columns[0])
    # Get the PLIF Tanimoto values in the row corresponding to the reference model
    plif_tanimoto = plif_tanimoto_df.loc[ref_model]
    plif_tanimoto = plif_tanimoto.tolist()
    # Add the PLIF Tanimoto values to the dataframe
    summary_df["plif_tanimoto"] = plif_tanimoto
    print(summary_df)

    # Prompt user for an output directory
    output_dir = input("Enter the output directory: ")
    # Create the output file path using the ligand name
    output_file = os.path.join(output_dir, ligand + "_summary.csv")
    # Write the dataframe to a csv file
    summary_df.to_csv(output_file, index=False)
    # Print a summary message to the user
    print("The summary file ({}) has been created in the specified directory.".format(os.path.basename(output_file)))

get_summary()
