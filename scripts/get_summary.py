"""
 This script will create a summary file containing the binding affinity, 
 PPS-Score, ligand RMSD, and PLIF Tanimoto values for each binding model. 
 The user will be prompted to enter the following information:                                                                           
                                                                                                  
 1. The directory containing the AutoDock Vina output files                                          
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
            "ensemble",
            "binding_affinity",
            "ppsscore",
            "lig_rmsd",
            "plif_tanimoto",
        ]
    )
    
    # Prompt user for directory containing the AutoDock Vina output files
    vina_logs = input(
        "Enter the path to the 'vina_output' directory: "
    )
    
    # Prompt user for ligand name
    ligand = input(
        "Enter the name of the ligand: "
    )
    
    # Initialize lists to hold binding affinities, species names, and ensemble numbers
    binding_affinities = []
    species = []
    ensembles = []
    file_count = 0

    # Find all *_ligand.pdbqt files in the vina_output directory
    # These files contain the docked ligand poses with binding affinity information
    vina_files = [
        f for f in os.listdir(vina_logs)
        if f.endswith("_ligand.pdbqt")
    ]
    
    # Sort files to ensure consistent ordering
    vina_files.sort()

    print(f"Found {len(vina_files)} ligand files to process")

    # Read binding affinity from each *_ligand.pdbqt file
    for file in vina_files:
        file_count += 1
        filepath = os.path.join(vina_logs, file)
        
        # Read the PDBQT file and extract binding affinity
        with open(filepath, "r") as f:
            lines = f.readlines()
            # Look for the binding affinity in the REMARK VINA RESULT line
            for line in lines:
                if line.startswith("REMARK VINA RESULT:"):
                    # Extract the binding affinity (4th element when split by whitespace)
                    binding_affinity = line.split()[3]
                    binding_affinities.append(binding_affinity)
                    
                    # Extract model name from filename (remove '_ligand.pdbqt' suffix)
                    model_name = file.replace("_ligand.pdbqt", "")
                    
                    # Parse species and ensemble from model name
                    # Split by underscore to get name components
                    name_parts = model_name.split("_")
                    
                    # Species is the first element
                    species_name = name_parts[0]
                    species.append(species_name)
                    
                    # Ensemble number is the last element if it's numeric, otherwise 0
                    # Examples:
                    #   "Chicken_AR_modified_ensemble_1" -> ensemble = 1
                    #   "Chicken_AR_modified" -> ensemble = 0
                    last_element = name_parts[-1]
                    if last_element.isdigit():
                        ensemble_num = int(last_element)
                    else:
                        ensemble_num = 0
                    ensembles.append(ensemble_num)
                    
                    # Only take the first REMARK VINA RESULT line (top pose)
                    break

    # Verify we got data for all files
    if len(binding_affinities) != len(vina_files):
        print(f"Warning: Expected {len(vina_files)} binding affinities but found {len(binding_affinities)}")

    # Add the binding affinities, species names, and ensemble numbers to the dataframe
    summary_df["binding_affinity"] = binding_affinities
    summary_df["species"] = species
    summary_df["ensemble"] = ensembles

    # Get the number of rows in the dataframe
    n_rows = len(summary_df.index)
    binding_models = list(range(1, n_rows + 1))
    # Add the binding models to the dataframe
    summary_df["binding_model"] = binding_models
    
    # Print the df
    print("\nBinding affinity data:")
    print(summary_df)

    # Prompt user for the directory containing the PPS-Score files
    ppsscore_dir = input(
        "\nEnter the path to the 'PPS_files' directory: "
    )
    
    # Create a list to hold the PPS-Score values
    ppsscores = []

    # Read the PPS-Score files
    # Sort to ensure consistent ordering with other data
    pps_files = sorted(os.listdir(ppsscore_dir))
    for file in pps_files:
        with open(os.path.join(ppsscore_dir, file), "r") as f:
            lines = f.readlines()
            # PPS-Score is on the third line (index 2)
            pps_line = lines[2]
            pps_line = pps_line.split()
            ppsscores.append(pps_line[2])

    # Add the PPS-Score values to the dataframe
    summary_df["ppsscore"] = ppsscores
    
    # Print the df
    print("\nAfter adding PPS-Scores:")
    print(summary_df)

    # Prompt user for the ligand RMSD file
    lig_rmsd_file = input("\nEnter the file path for the ligand RMSD file: ")
    # Remove quotes from the file path
    lig_rmsd_file = lig_rmsd_file.replace('"', "")
    
    # Read the ligand RMSD file
    with open(lig_rmsd_file, "r") as f:
        lines = f.readlines()
        lines = [line.split() for line in lines]
        # Skip the first two header lines
        lig_rmsds = [line[1] for line in lines[2:]]
    
    # Add the ligand RMSD values to the dataframe
    summary_df["lig_rmsd"] = lig_rmsds
    
    # Print the df
    print("\nAfter adding ligand RMSDs:")
    print(summary_df)

    # Prompt user for the PLIF Tanimoto file
    plif_tanimoto_file = input(
        "\nEnter the file path for the plif_similarity_summary.csv file: "
    )
    plif_tanimoto_file = plif_tanimoto_file.replace('"', "")
    plif_tanimoto_df = pd.read_csv(plif_tanimoto_file)
    
    # Sort the PLIF Tanimoto dataframe by test_structure to ensure consistent ordering
    plif_tanimoto_df = plif_tanimoto_df.sort_values(by="test_structure")
    
    # Copy the PLIF Tanimoto values to the summary dataframe
    summary_df["plif_tanimoto"] = plif_tanimoto_df["tanimoto_similarity"].values
    
    print("\nFinal summary with all data:")
    print(summary_df)

    # Prompt user for an output directory
    output_dir = input("\nSpecify the output directory: ")
    
    # Create the output file path using the ligand name
    output_file = os.path.join(output_dir, ligand + "_summary.csv")
    
    # Write the dataframe to a csv file
    summary_df.to_csv(output_file, index=False)
    
    # Print a summary message to the user
    print(
        "\nThe summary file ({}) has been created in the specified directory.".format(
            os.path.basename(output_file)
        )
    )


if __name__ == "__main__":
    get_summary()