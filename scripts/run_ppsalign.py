import os
import subprocess


def ppsalign_loop():
    """
    This function will take a file directory of .poc files as an input and 
    run PPSalign on each.     
    The user will need to ensure that the PPSalign executable is in their PATH.

    Input: .poc files
    Output: .txt files containing PPS-Score values
    """
    # Get the path to the directory containing the .poc files
    poc_dir = input("Enter the path to the 'binding_sites' directory containing the .poc files: ")
    # Check if the path exists and if not, ask for it again
    while not os.path.exists(poc_dir):
        poc_dir = input(
            "That path does not appear to exist.\nPlease enter the path to the directory containing the .poc files: "
        )
    # Change the working directory to the directory containing the .poc files
    os.chdir(poc_dir)
    # Get the list of .poc files in the directory
    poc_files = [i for i in os.listdir(poc_dir) if i.endswith(".poc")]

    if not poc_files:
        print("No .poc files found in the directory.")
        return

    # Search for the template .poc file
    template = None
    for poc in poc_files:
        if poc.startswith("ref_"):
            ref = input(f"Is {poc} the template .poc file? (y/n): ").lower()
            while ref not in ["y", "n"]:
                ref = input("Please enter y or n: ").lower()
            if ref == "y":
                template = poc
                break

    if not template:
        template = input("Enter the name of the template (reference) .poc file: ")
        if not template.endswith(".poc"):
            template += ".poc"
        while template not in poc_files:
            template = input(
                "That file does not appear to exist.\nPlease enter the name of the template (reference) .poc file: "
            )
            if not template.endswith(".poc"):
                template += ".poc"

    # Create output directory if it doesn't exist
    if not os.path.isdir("PPS_files"):
        os.mkdir("PPS_files")

    for poc in poc_files:
        # Skip the template file
        if poc != template:
            output_file = os.path.join("PPS_files", poc.split(".poc")[0] + "_PPS.txt")
            print("Calculating PPS-Score for " + poc + "...")
            # Run PPSalign with proper argument list
            try:
                # Run PPSalign without check=True
                subprocess.run(
                    ["PPSalign", poc, template],
                    stdout=open(output_file, "w")
                )
            except FileNotFoundError:
                print("Error: PPSalign command not found. Make sure it's in your PATH.")
                return

    print("PPS-Score calculation completed and files were saved to 'PPS_files' folder.")


if __name__ == "__main__":
    ppsalign_loop()