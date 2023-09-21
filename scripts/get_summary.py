"""
 This script will create a summary file containing the binding affinity, 
 PPS-Score, ligand RMSD, and PLIF Tanimoto values for each binding model. 
 The user will be prompted to enter the following information:                                                                           
                                                                                                  
 1. The directory containing the AutoDock Vina log files                                          
 2. The name of the ligand, which should match the suffix of the log files                        
 3. The file path for the PPS-Score file                                                          
 4. The file path for the ligand RMSD file                                                        
 5. The file path for the PLIF Tanimoto matrix file                                               
 6. The name of the reference model (as listed in the matrix)                                     
 7. The output directory                                                                          
                                                                                                  
 The summary file will be written to the specified output directory.                              
"""

import os

import pandas as pd


def get_summary():
    """
    Main function.
    Prompts user for the necessary information and creates the summary file.
    """
    # Initialize empty dataframe to store the summary data
    summary_df = pd.DataFrame(
        columns=[
            "binding_model",
            "species",
            "pose",
            "binding_affinity",
            "ppsscore",
            "lig_rmsd",
            "plif_tanimoto",
        ]
    )
    # Prompt user for directory containing the AutoDock Vina log files
    vina_logs = input(
        "Enter the file directory containing the (split) AutoDock Vina log files: "
    )
    # Prompt user for ligand name
    ligand = input(
        "Enter the name of the ligand: "
    )
    # Initialize lists to hold binding affinities, poses, and species names
    binding_affinities = []
    poses = []
    species = []
    file_count = 0

    vina_dir = [
        i
        for i in os.listdir(vina_logs)
        if i.endswith(".pdbqt") and "bound" in i and "ligand" in i
    ]

    # Read the AutoDock Vina log files
    for file in vina_dir:
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
                    # Get the species name by splitting the file name
                    # at the underscore
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

    # Ask user if there is a combined PPS-Score file or a directory of PPS-Score files
    ppsscore_type = input(
        "Is there a combined PPS-Score file or a directory of PPS-Score files? "
        "(Enter 'combined' or 'directory'): "
    )

    while ppsscore_type not in ["combined", "directory"]:
        ppsscore_type = input("Please enter 'combined' or 'directory': ")

    # If there is a combined PPS-Score file
    if ppsscore_type == "combined":
        # Prompt user for the path to the combined PPS-Score file
        ppsscore_file = input("Enter the file path for the combined PPS-Score file: ")
        ppsscore_file = ppsscore_file.replace('"', "")
        # Create a list to hold the PPS-Score values
        ppsscores = []

        # Read the PPS-Score file
        with open(ppsscore_file, "r") as f:
            lines = f.readlines()
            line_counter = 2
            for line in lines[2 : file_count + 2]:
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

    # If there is a directory of PPS-Score files
    elif ppsscore_type == "directory":
        # Prompt user for the directory containing the PPS-Score files
        ppsscore_dir = input(
            "Enter the file directory containing the PPS-Score files: "
        )
        # Create a list to hold the PPS-Score values
        ppsscores = []

        # Read the PPS-Score files
        for file in os.listdir(ppsscore_dir):
            with open(os.path.join(ppsscore_dir, file), "r") as f:
                lines = f.readlines()
                pps_line = lines[2]
                pps_line = pps_line.split()
                ppsscores.append(pps_line[2])

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
    plif_tanimoto_file = input(
        "Enter the file path for the PLIF Tanimoto matrix .csv file: "
    )
    plif_tanimoto_file = plif_tanimoto_file.replace('"', "")
    plif_tanimoto_df = pd.read_csv(plif_tanimoto_file)

    # Search for the reference model in the PLIF Tanimoto matrix file
    ref_model = ""
    # Loop through each column name in the dataframe
    for col in plif_tanimoto_df.columns:
        # If the column name starts with "ref_", ask the user if it is the reference model
        if col.startswith("ref_"):
            ref_model = input(
                "Is {} the reference model? (y/n): ".format(col)
            )
            ref_model = ref_model.lower()
            # Check if the user entered a valid response
            while ref_model not in ["y", "n"]:
                ref_model = input(
                    "Please enter 'y' or 'n': "
                )
                ref_model = ref_model.lower()
            # If the user says yes, break out of the loop
            if ref_model == "y":
                ref_model = col
                break
    if ref_model == "n":
        ref_model = input(
            "Enter the name of the reference model (as listed in the matrix): "
        )

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
    output_dir = input("Specify the output directory: ")
    # Create the output file path using the ligand name
    output_file = os.path.join(output_dir, ligand + "_summary.csv")
    # Write the dataframe to a csv file
    summary_df.to_csv(output_file, index=False)
    # Print a summary message to the user
    print(
        "The summary file ({}) has been created in the specified directory.".format(
            os.path.basename(output_file)
        )
    )


if __name__ == "__main__":
    get_summary()
