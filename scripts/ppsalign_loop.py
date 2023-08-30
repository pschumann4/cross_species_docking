"""
 This function will take a file directory of .poc files as an input and 
 run PPSalign on each.     
 The user will need to set "PPSalign" as an environment variable.                                 
"""

import os


def ppsalign_loop():
    """
    Main funciton
    """
    # Get the path to the directory containing the .poc files
    poc_dir = input("Enter the path to the directory containing the .poc files: ")
    # Check if the path exists and if not, ask for it again
    while not os.path.exists(poc_dir):
        poc_dir = input(
            "That path does not appear to exist.\nPlease enter the path to the directory containing the .poc files: "
        )
    # Change the working directory to the directory containing the .poc files
    os.chdir(poc_dir)
    # Get the list of .poc files in the directory
    poc_files = [i for i in os.listdir(poc_dir) if i.endswith(".poc")]

    # Search for the template .poc file
    template = ""
    for poc in poc_files:
        if poc.startswith("ref_"):
            ref = input("Is {} the template .poc file? (y/n): ".format(poc))
            ref = ref.lower()
            # check for valid input
            while ref.lower() not in ["y", "n"]:
                ref = input("Please enter y or n: ")
            if ref == "y":
                template = poc
                break
    if ref == "n":
        template = input("Enter the name of the template (reference) .poc file: ")
        if not template.endswith(".poc"):
            template += ".poc"
        while template not in poc_files:
            template = input(
                "That file does not appear to exist.\nPlease enter the name of the template (reference) .poc file: "
            )
            if not template.endswith(".poc"):
                template += ".poc"

    for poc in poc_files:
        if poc != template:
            # Create the command to run PPSalign
            cmd = (
                "PPSalign "
                + poc
                + " "
                + template
                + " > "
                + poc.split(".poc")[0]
                + "_PPS.txt"
            )
            print("Calculating PPS-Score for " + poc + "...")
            os.system(cmd)
        else:
            continue

    # Move all of the _PPS.txt files to a folder called 'PPS_files'
    if not os.path.isdir("PPS_files"):
        os.mkdir("PPS_files")
    for file in os.listdir(poc_dir):
        if file.endswith("_PPS.txt"):
            os.rename(file, "PPS_files/" + file)
    print("Done!")


if __name__ == "__main__":
    ppsalign_loop()
