import os
import subprocess
import glob
import shutil


def run_vina():
    """
    Run AutoDock Vina for each configuration file in the specified directory.
    Returns True if successful, False otherwise.
    """
    try:
        # Find all configuration files
        config_files = glob.glob("*_conf.txt")
        
        if not config_files:
            print("No configuration files found in the specified directory.")
            return False
        
        # Process each configuration file
        for config_file in config_files:
            try:
                # Get base name (removing _conf.txt)
                base_name = config_file.replace("_conf.txt", "")
                
                print(f"\nRunning AutoDock Vina on '{base_name}'...")
                
                # Run Vina with the configuration file
                result = subprocess.run(
                    ["vina", "--config", config_file],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                # Print Vina output
                print(result.stdout)
                
                if result.stderr:
                    print("Errors/Warnings:")
                    print(result.stderr)
                    
            except subprocess.CalledProcessError as e:
                print(f"Error running Vina for {config_file}:")
                print(e.stderr)
                return False
            except Exception as e:
                print(f"Unexpected error processing {config_file}:")
                print(str(e))
                return False
        
        return True
        
    except Exception as e:
        print(f"Unexpected error in run_vina: {str(e)}")
        return False

def process_vina_output(chemical):
    """
    Process Vina output files by splitting them and organizing the results.
    Creates an output directory and moves split files there.
    Returns True if successful, False otherwise.
    """
    print("\nRunning vina_split on Vina output files...")
    
    try:
        # Create output folder in current directory
        output_folder = "vina_output"
        os.makedirs(output_folder, exist_ok=True)
        
        # Find all relevant PDBQT files
        pattern = f"*bound_{chemical}"
        pdbqt_files = glob.glob(pattern)
        
        if not pdbqt_files:
            print(f"No files matching pattern '{pattern}' found in the specified directory.")
            return False
        
        # Process each file
        for pdbqt_file in pdbqt_files:
            try:
                print(f"Processing file: {pdbqt_file}")
                
                # Run vina_split
                result = subprocess.run(
                    ["vina_split", "--input", pdbqt_file],
                    capture_output=True,
                    text=True,
                    check=True
                )
                
                if result.stderr:
                    print("Warnings/Errors from vina_split:")
                    print(result.stderr)
                
                # Get base name without extension
                base_name = os.path.splitext(pdbqt_file)[0]
                
                # Move split files to output directory
                for i in range(1, 101):  # Range 1-100
                    patterns = [
                        f"{base_name}_flex_{i}.pdbqt",
                        f"{base_name}_rigid_{i}.pdbqt",
                        f"{base_name}_ligand_{i}.pdbqt"
                    ]
                    
                    # Try to move each possible file
                    for pattern in patterns:
                        if os.path.exists(pattern):
                            try:
                                shutil.move(pattern, os.path.join(output_folder, pattern))
                            except Exception as e:
                                print(f"Error moving {pattern}: {str(e)}")
                                return False
                
            except subprocess.CalledProcessError as e:
                print(f"Error running vina_split for {pdbqt_file}:")
                print(e.stderr)
                return False
            except Exception as e:
                print(f"Unexpected error processing {pdbqt_file}:")
                print(str(e))
                return False
        
        print("Processing completed successfully.")
        return True
        
    except Exception as e:
        print(f"Unexpected error in process_vina_output: {str(e)}")
        return False

def run_vina_batch():
    """
    Run AutoDock Vina for a batch of PDBQT files in a specified directory.
    This requires that the specied directory contains:
    1. The protein receptor file in PDBQT format.
    2. The ligand file in PDBQT format.
    3. Configuration files for each protein-ligand pair.

    This function outputs the Vina results in a new directory called 'vina_output'.
    """
    # Get the configuration directory from user
    pdbqt_dir = input("Enter the path to the PDBQT file directory: ")

    # Get the name of the chemical file used for the binding simulation
    chemical = input("What is the name of the chemical file used for the binding simulation: ")
    if not chemical.endswith(".pdbqt"):
        chemical += ".pdbqt"
    while not chemical:
        chemical = input("That file does not appear to exist.\n"
                            "Please try entering the name again (case sensitive): ")
        if not chemical.endswith(".pdbqt"):
            chemical += ".pdbqt"
    
    # Verify the directory exists
    while not os.path.exists(pdbqt_dir):
        pdbqt_dir = input("That path does not appear to exist.\n"
                          "Please enter the path to the PDBQT file directory: ")

    # Change to the configuration directory
    os.chdir(pdbqt_dir)
    
    # Run Vina
    if run_vina():
        # Process Vina output
        process_vina_output(chemical)

if __name__ == "__main__":
    run_vina_batch()