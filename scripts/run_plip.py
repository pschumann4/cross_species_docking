import os
import shutil
import subprocess


def process_pdb_files():
    """
    Run PLIP on all PDB files in a specified directory and move the results to a new folder.
    OpenBabel is used to process protonated PDB files by adding hydrogens for consistency.
    """
    # Get directory path from user
    dir_path = input("Enter the path to the folder with the PDB files to analyze: ")
    
    # Change to specified directory
    os.chdir(dir_path)
    
    # Run PLIP on all PDB files
    for file in os.listdir(dir_path):
        if file.endswith(".pdb"):
            base_name = os.path.splitext(file)[0]
            try:
                print(f"Running PLIP on {file}...")
                subprocess.run(["plip", "-f", file, "-xv", "--name", base_name])
            except Exception as e:
                print(f"Error processing {file}: {e}")
    
    # Process protonated files with OpenBabel
    for file in os.listdir(dir_path):
        if file.endswith('_protonated.pdb'):
            subprocess.run(["obabel", file, "-o", "pdb", "-O", file, "-h"])
    
    # Create and move files to plip_results folder
    results_dir = os.path.join(dir_path, "plip_results")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    
    print("Moving files to the plip_results folder...")
    
    # Move XML and protonated PDB files
    for file in os.listdir(dir_path):
        if file.endswith(".xml") or "_protonated.pdb" in file:
            src = os.path.join(dir_path, file)
            dst = os.path.join(results_dir, file)
            shutil.move(src, dst)

if __name__ == "__main__":
    process_pdb_files()