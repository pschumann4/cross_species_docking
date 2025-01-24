import os
import subprocess


def prepare_receptor(pdb, pdb_name, flexible_residues=None):
    """
    Prepare the receptor using mk_prepare_receptor command
    
    Parameters:
    -----------
    pdb : str
        Path to the input PDB file
    pdb_name : str
        Name of the PDB file without extension
    flexible_residues : str, optional
        Flexible residues to prepare, formatted as A:1,A:2,...
    """
    command = ['mk_prepare_receptor', '--read_pdb', pdb, '-o', pdb_name, '-p', '-a']
    if flexible_residues:
        command.extend(['-f', flexible_residues])
    
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError:
        print(f"\nFailed to prepare the receptor for {pdb_name}.")
        if flexible_residues:
            print("\nConsider preparing the receptor using the PDBFixer and try again.")
    except FileNotFoundError:
        print("\nThe command 'mk_prepare_receptor' was not found. Please ensure it is installed and available in your PATH.")

def get_flexible_residues(residues_file):
    """
    Read flexible residues from a file
    
    Parameters:
    -----------
    residues_file : str
        Path to the flex_residues.txt file
    
    Returns:
    --------
    dict: Dictionary with PDB names as keys and flexible residues as values
    """
    with open(residues_file, "r") as f:
        res_list = f.readlines()[1:]  # Remove header line
        res_list = [line for line in res_list if line.strip()]  # Remove empty lines
    
    flex_residues = {}
    for line in res_list:
        pdb_name, res = line.split(":")
        res = res.split("[")[1].split("]")[0]
        flex_residues[pdb_name] = res
    
    return flex_residues

def main():
    pdb_dir = input("Enter the path to the directory containing the aligned PDB files: ")
    while not os.path.exists(pdb_dir):
        print("That path does not appear to exist.")
        pdb_dir = input("\nPlease enter the path to the directory containing the aligned PDB files: ")
    
    os.chdir(pdb_dir)
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(".pdb")]

    flex_residues = input("Are there flexible residues to prepare? (y/n): ").lower()
    while flex_residues not in ["y", "n"]:
        flex_residues = input("Please enter y or n: ").lower()

    if flex_residues == "y":
        residues_file = "flex_residues.txt"
        if not os.path.exists(residues_file):
            residues_file = input("Enter the path to the flex_residues.txt file: ").strip('"')
            while not os.path.exists(residues_file):
                residues_file = input("That path does not appear to exist.\nPlease enter the path to the flex_residues.txt file: ").strip('"')
        
        flex_residues_dict = get_flexible_residues(residues_file)
    else:
        flex_residues_dict = {}

    for pdb in pdb_files:
        if pdb.startswith("ref_"):
            print(f"\nSkipping reference structure {pdb}...")
            continue
        
        pdb_name = pdb.split(".pdb")[0]
        print(f"\nPreparing {pdb_name}...")

        if pdb_name in flex_residues_dict:
            res = flex_residues_dict[pdb_name]
            if res:
                res_id = [i[3:] for i in res.split(",")]
                with open(pdb, "r") as f:
                    pdb_lines = f.readlines()
                    for line in pdb_lines:
                        if line[22:26].strip() in res_id:
                            chain = line[21]
                            res_id[res_id.index(line[22:26].strip())] = chain + ":" + line[22:26].strip()
                res_id = ",".join(res_id)
                prepare_receptor(pdb, pdb_name, res_id)
            else:
                print(f"No flexible residues for {pdb_name}. Preparing rigid receptor...")
                prepare_receptor(pdb, pdb_name)
        else:
            prepare_receptor(pdb, pdb_name)

    if not os.path.exists("pdbqt_files"):
        os.mkdir("pdbqt_files")
    for f in os.listdir(pdb_dir):
        if f.endswith(".pdbqt"):
            os.replace(f, os.path.join("pdbqt_files", f))
    print("\nSuccessfully prepared structures have been moved to the 'pdbqt_files' directory.")

if __name__ == "__main__":
    main()