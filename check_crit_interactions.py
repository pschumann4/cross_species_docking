####################################################################################################
# This script allows the user to check for specific interactions involving critical residues in a 
# set of PLIF text files. The user defined residues will be highlighted in a heatmap.
####################################################################################################

############################
# Import required packages #
############################

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

######################
# Load all functions #
######################

def check_exit(input_str):
    if isinstance(input_str, str):
        input_str = input_str.lower()
        if input_str == 'exit':
            return True
        return False

def check_crit_interactions():
    # Specify a dictionary of valid residue names including genetically encoded non-cannonical amino acids
    aa_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'SEC': 'U', 'PYL': 'O'}    # Create lists to store the residue info
    restype = []
    resnr = []

    # Prompt user for file directory containing the PLIF text files
    pwd = input("Enter the file directory of the PLIF text file: ")
    if check_exit(pwd):
        return
    
    # Check that the directory exists
    while not os.path.exists(pwd):
        print('Error: This directory does not exist. Make sure that the directory path is correct.\n')
        pwd = input('Enter your working directory: ')
        if check_exit(pwd):
            return
    # Change the working directory to the specified directory
    os.chdir(pwd)

    # Prompt user for comma separated list of residue names and numbers
    crit_res = input("Enter a comma separated list of residue names and numbers (e.g. 'ARG 1, PHE 2, LEU 3'): ")
    if check_exit(crit_res):
        return
    # Split the string into a list of residue names and numbers
    if ", " in crit_res:
        crit_res = crit_res.split(", ")
    elif "," in crit_res:
        crit_res = crit_res.split(",")
    # Capitalize the residue name
    for i in range(len(crit_res)):
        crit_res[i] = crit_res[i].upper()

    # Check that a space was entered between the residue name and number
    for i in range(len(crit_res)):
        if " " not in crit_res[i]:
            print("Please enter a space between the residue name and number.")
            crit_res = input("Enter a comma separated list of residue names and numbers (e.g. 'ARG 1, PHE 2, LEU 3'): ")
            if check_exit(crit_res):
                return
            if ", " in crit_res:
                crit_res = crit_res.split(", ")
            elif "," in crit_res:
                crit_res = crit_res.split(",")
            for i in range(len(crit_res)):
                crit_res[i] = crit_res[i].upper()

    # Split the residue names and numbers into separate lists
    for i in range(len(crit_res)):
        crit_res[i] = crit_res[i].split(" ")
        restype.append(crit_res[i][0])
        resnr.append(crit_res[i][1])

    # If the user entered single letter amino acid codes, convert them to three letter codes
    for i in range(len(restype)):
        if len(restype[i]) == 1:
            for key, value in aa_dict.items():
                if restype[i] == value:
                    restype[i] = key

    # Check that all of the residue names are valid
    for i in range(len(restype)):
        if restype[i] not in aa_dict.keys():
            print("Invalid residue name: " + restype[i] + "\nPlease enter a valid residue name.")
            return
    # Create a dictionary to store the results
    results = {}

    # Loop through the PLIF files in the directory
    for file in os.listdir(pwd):
        if file.endswith(".txt"):
            # Initialize list to store the temporary results
            conserved_res = []
            # Get the base name of the file
            base = os.path.basename(file).split("_PLIF")[0]
            # Read the PLIF text file and convert it to a pandas dataframe
            df = pd.read_csv(file, sep = "\s+")
            # Search each row of the dataframe for the critical residue(s)
            # If the 'query' column of that residue is '1', add the residue to the list
            for i in range(len(df)):
                # Loop over each critical residue
                for j in range(len(restype)):
                    # If the residue name and number match, and the 'query' column is '1', add the residue to the list
                    if df['restype'][i] == restype[j] and df['resnr'][i] == int(resnr[j]) and df['query'][i] == 1:
                        # Check that the residue is not already in the list
                        if restype[j] + " " + resnr[j] not in conserved_res:
                            conserved_res.append(restype[j] + " " + resnr[j])
            # Add the results to the dictionary with the file name as the key
            results[base] = conserved_res

    # Create a pandas dataframe from the dictionary
    df = pd.DataFrame.from_dict(results, orient = 'index')
    df2 = pd.DataFrame(index=df.index)

    for i, row in df.iterrows():
    # Loop through each element of the row, starting from the index 0
        for j in range(0, len(row)):
            # If the element is not None, split it into amino acid and position
            if row[j] is not None:
                aa, pos = row[j].split()
                # Create a new column with the name of the amino acid and position
                col_name = f'{aa} {pos}'
                # Set the value of the cell in the new dataframe to 1
                df2.loc[i, col_name] = 1

    # Fill the NaN values with 0
    df2 = df2.fillna(0)
    # Write the dataframe to a csv file
    df2.to_csv("crit_res_interactions.csv")
    # Convert the dataframe to a heatmap using a green-red color scheme
    sns.heatmap(df2, cmap = 'RdYlGn', annot = False, cbar = False, linecolor='black', linewidths=1)
    plt.tight_layout()  
    # Save the heatmap as a png file
    plt.savefig("crit_res_interactions.png", dpi = 300)
    # Display the heatmap
    plt.show()
    return

check_crit_interactions()