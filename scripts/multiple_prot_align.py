"""
 This script is used to update the residue numbers of multiple proteins using a               
 multiple sequence alignment (MUSCLE) and to update the protein coordinates using 
 structural alignments via PyMOL. Additionally, this function will remove extraneous 
 chains from the PDB files and remove duplicate atoms, if requested. The aligned 
 structures will also be trimmed to include only the residues that are aligned 
 to the reference structure (+10 residue buffer at each terminus).                                                       
 The script will output a text file containing the RMSD values for each structure
 as well as a csv file containing the original residue numbers and the new residue 
 numbers. Lastly, the aligned and modified structures will be saved as PDB files 
 along with copies of the original PDB files.                                                                   
                                                                                              
 The user will need to have PyMOL installed on their computer and have the MUSCLE 
 executable added to their PATH as "muscle".                                                             
"""

import os
import sys
import shutil
import subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import check_tools

import pandas as pd
from Bio import AlignIO

# Check that the user has PyMOL installed on their computer
try:
    import pymol
    from pymol import cmd
except ImportError:
    print(
        "PyMOL is not installed on this computer. Please install PyMOL before running this script."
    )


def pdb_to_fasta(pdb_file):
    """
    This function will take a PDB file and return the FASTA sequence for the protein.
    """
    # Create a dictionary to convert the three letter amino acid codes to one letter codes
    aa_dict = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLN": "Q",
        "GLU": "E",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "SEC": "U",
        "PYL": "O",
    }
    # Create a list to store the FASTA sequence
    fasta_seq = []
    # Create a list to track the residue numbers to avoid duplicates
    resnum_list = []
    # Open the PDB file
    with open(pdb_file, "r") as pdb:
        # Iterate over the lines in the PDB file
        for line in pdb:
            # Check if the line is a standard amino acid residue
            if line.startswith("ATOM") and line[17:20] in aa_dict.keys():
                # Get the residue number and the residue name
                resnum = line[22:27].strip()
                resname = line[17:20]
                # Check if the residue number is already in the list
                if resnum not in resnum_list:
                    # Add the residue number to the list
                    resnum_list.append(resnum)
                    # Add the one letter code to the FASTA sequence
                    fasta_seq.append(aa_dict[resname])
    # Join the FASTA sequence into a string
    fasta_seq = "".join(fasta_seq)
    return fasta_seq


def pdb_dir_to_fasta_dict(pwd):
    """
    Converts all PDB files in a directory to FASTA sequences and returns a dictionary
    of the FASTA sequences.
    """
    # Create an empty dictionary to store the FASTA sequences
    fasta_dict = {}
    files = [i for i in os.listdir(pwd) if i.endswith(".pdb")]
    # Loop through the files in the directory
    for file in files:
        # Get the name of the PDB file without the extension
        pdb_name = os.path.splitext(file)[0]
        # Get the FASTA sequence for the PDB file
        fasta_seq = pdb_to_fasta(os.path.join(pwd, file))
        # Add the FASTA sequence to the dictionary
        fasta_dict[pdb_name] = fasta_seq
    return fasta_dict


def align_sequences(fasta_dict, output_file):
    """
    Performs a multiple sequence alignment on a dictionary of FASTA sequences and
    returns a dictionary of the aligned sequences.
    """
    # Write the FASTA sequences to a file using their names as the
    # sequence identifiers
    with open("input.fasta", "w") as fasta_file:
        for key, value in fasta_dict.items():
            fasta_file.write(">" + key + "\n" + value + "\n")
    # Run MUSCLE on the input.fasta file
    subprocess.run(
        ["muscle", "-align", "input.fasta", "-output", output_file],
        check=True,
        capture_output=True
    )
    # Read the output file
    aligned = AlignIO.read(output_file, "fasta")
    # Create an empty dictionary to store the aligned sequences
    aligned_dict = {}
    for record in aligned:
        aligned_dict[record.id] = str(record.seq)
    return aligned_dict


def get_aligned_residues(aligned_dict):
    """
    Generates a dictionary of the aligned residues for each sequence in the alignment.
    """
    # Create an empty dictionary to store the new residue positions for each sequence
    new_positions = {}
    # Start a counter at 1 to track the residue position
    position = 1
    # Loop through the sequences in the alignment dictionary
    for key, value in aligned_dict.items():
        # Create a list to store the residue positions
        position_list = []
        # Loop through the residues in the sequence
        for residue in value:
            # Check if the residue is a gap, and if so, increment the position counter
            if residue == "-":
                position += 1
            # If the residue is not a gap, add the position to the list
            # and increment the position counter
            else:
                position_list.append(position)
                position += 1
        # Add the list of positions to the dictionary
        new_positions[key] = position_list
        # Reset the position counter
        position = 1
    return new_positions


def modify_pdb_sequence(pdb_file, positions, output_file):
    """
    Modify the sequence of a PDB file to match the aligned sequence.
    """
    # Create an empty list to store the modified records
    atom_records = []
    hetatom_records = []
    ter_records = []
    other_records = []
    with open(pdb_file, "r") as f:
        # Initialize a counter to keep track of the current position
        position_index = 0
        # Create a dictionary to store the residue names and numbers
        residues = {}
        for line in f:
            if line.startswith("ATOM"):
                # Extract the relevant information from the ATOM record
                atom_number = line[6:11]
                atom_name = str(line[12:16]).rjust(4)
                residue_name = line[17:20]
                chain_id = line[21]
                residue_number = line[22:26]
                x = line[30:38]
                y = line[38:46]
                z = line[46:54]
                occupancy = line[54:60]
                temperature_factor = line[60:66]
                element = line[76:78]
                # Check if the residue has been encountered before
                if residue_number not in residues:
                    # If the residue is new, add it to the dictionary
                    residues[residue_number] = positions[position_index]
                    position_index += 1
                # Modify the residue number
                new_residue_number = str(residues[residue_number]).rjust(4)
                # Create a new ATOM record with the modified residue number and position
                modified_record = f"ATOM  {atom_number} {atom_name} {residue_name} {chain_id}{new_residue_number}    {x}{y}{z}{occupancy}{temperature_factor}          {element}\n"
                atom_records.append(modified_record)
            else:
                if line.startswith("HETATM"):
                    hetatom_records.append(line)
                if line.startswith("TER"):
                    # If the TER record is empty, add it to the list. Otherwise, modify it.
                    if line.strip() == "TER":
                        ter_records.append(line)
                        continue
                    else:
                        record_type = line[0:6]
                        serial_number = line[6:11]
                        residue_name = line[17:20]
                        chain_id = line[21]
                        # Replace the residue number with the last element of the 'positions' list
                        new_residue_number = str(positions[-1]).rjust(4)
                        # Create a new TER record with the modified residue number
                        modified_ter = f"{record_type}{serial_number}      {residue_name} {chain_id}{new_residue_number}\n"
                        ter_records.append(modified_ter)
                if line.startswith("END"):
                    ter_records.append(line)
                if (
                    line.startswith("CRYST")
                    or line.startswith("SCALE")
                    or line.startswith("ORIGX")
                    or line.startswith("REMARK")
                    or line.startswith("MODEL")
                ):
                    other_records.append(line)
    # Remove duplicate ter_records
    ter_records = list(dict.fromkeys(ter_records))
    # Reorder the modified 'ATOM' records by the integer value of the residue number
    atom_records = sorted(atom_records, key=lambda x: int(x[22:27]))
    # Combine the modified records into a single list
    modified_records = other_records + hetatom_records + atom_records + ter_records
    # Create the output file directory as the same as the input file directory
    output_file = os.path.join(os.path.dirname(pdb_file), output_file)
    # Write the modified records to the output file directory
    with open(output_file, "w") as f:
        for record in modified_records:
            f.write(record)
    return output_file


def align_structures(pwd, ref):
    """
    Performs as structural alignment of all pdb files in the working
    directory to the reference structure using PyMOL.
    """
    # Loop through the files in the working directory to get the pdb files that
    # are not the reference structure
    structures = []
    for file in os.listdir():
        if (
            file.endswith("modified.pdb")
            and file != ref
            and file != f"{ref[:-4]}_modified.pdb"
        ):
            structures.append(file)
    # Check that there are pdb files in the working directory that are not the
    # reference structure
    if len(structures) == 0:
        print("There are no other pdb files in the working directory.")
        return
    # Initialize PyMOL
    pymol.finish_launching()
    # Load the reference structure
    cmd.load(ref, "ref")
    # Remove waters
    cmd.remove("resn HOH")
    # Remove hydrogens
    cmd.remove("hydro")
    # Save the reference structure as a PDB file using the original name followed by "_ref"
    cmd.save(f"ref_{ref[:-4]}.pdb", "ref", 1)
    # Remove heteroatoms
    cmd.remove("hetatm")
    # Save the reference structure as a PDB file using the original name
    cmd.save(ref, "ref", 1)
    # Create a dictionary to store the RMSD values
    rmsd_dict = {}
    for structure in structures:
        # Load the structure
        cmd.load(structure, "s")
        # Remove all heteroatoms
        cmd.remove("hetatm")
        # Remove hydrogens
        cmd.remove("hydro")
        # Align the structure to the reference
        cmd.align("s", "ref")
        # Store the RMSD value
        rmsd = cmd.align("s", "ref")[0]
        rmsd = round(rmsd, 3)
        rmsd_dict[structure] = rmsd
        # save the aligned structure with "_aligned" appended to the name
        cmd.save(structure, "s", 1)
        cmd.delete("s")
    # View the aligned structures
    cmd.color("red", "ref")
    for structure in structures:
        # load the structure and use the structure name as the object name
        cmd.load(structure, structure)
    # Write the RMSD values to a file
    with open("prot_rmsd_values.txt", "w") as f:
        f.write("Reference structure: %s" % ref + "\n")
        for structure, rmsd in rmsd_dict.items():
            f.write("%s: %s\n" % (structure, rmsd))


def get_og_pos(og_pdb, new_pdb):
    """
    This function will get the original residue positions from the original pdb file.
    """
    # Get the base name of the first pdb file
    og_pdb_name = os.path.splitext(og_pdb)[0]
    # Read in the pdb files
    with open(og_pdb) as f1, open(new_pdb) as f2:
        og_pdb = f1.readlines()
        new_pdb = f2.readlines()
    # Create a list to store the original positions
    og_pos = []
    og_res = []
    # Create a list to store the new positions
    new_pos = []
    new_res = []
    for line in og_pdb:
        if line.startswith("ATOM"):
            res_name = line[17:20].strip()
            res_num = line[22:27].strip()
            if res_num in og_pos:
                continue
            og_pos.append(res_num)
            og_res.append(res_name)
    for line in new_pdb:
        if line.startswith("ATOM"):
            res_name = line[17:20].strip()
            res_num = line[22:27].strip()
            if res_num in new_pos:
                continue
            new_pos.append(res_num)
            new_res.append(res_name)
    # Create a dataframe to store the data
    df = pd.DataFrame(
        {"new_pos": new_pos, og_pdb_name + "_res": og_res, og_pdb_name + "_pos": og_pos}
    )
    return df


def remove_chains(pdb, keep_chain=None):
    """
    Removes extraneous chains from a PDB file.
    
    Args:
        pdb (str): Path to the PDB file
        keep_chain (str, optional): Specific chain ID to keep. If None, user will be prompted.
        
    Returns:
        str or None: Returns kept chain IDs if same_chain prompt is answered 'y', None otherwise
        
    Raises:
        ValueError: If specified keep_chain is not found in the PDB file
    """
    dir = os.path.dirname(pdb)
    filename = os.path.basename(pdb)
    
    def get_available_chains(lines):
        """Helper function to get all unique chain IDs from PDB file."""
        chains = set()
        for line in lines:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and len(line) >= 22:
                chains.add(line[21])
        return sorted(list(chains))
    
    def write_filtered_pdb(lines, keep_chains):
        """Helper function to write filtered PDB file with specified chains."""
        new_pdb = []
        for line in lines:
            # Keep non-coordinate records
            if not any(line.startswith(record) for record in ["ATOM", "HETATM", "TER", "ANISOU"]):
                if not line.startswith("END"):
                    new_pdb.append(line)
                continue
            
            # Process coordinate records
            if len(line) >= 22:
                chain = line[21]
                if isinstance(keep_chains, (list, set)) and chain in keep_chains:
                    new_pdb.append(line)
                elif isinstance(keep_chains, str) and chain == keep_chains:
                    new_pdb.append(line)
        
        # Ensure END record is present
        if new_pdb and not new_pdb[-1].startswith("END"):
            new_pdb.append("END\n")
            
        # Write the filtered PDB file
        with open(os.path.join(dir, filename), "w") as f:
            f.writelines(new_pdb)
        print(f"Extraneous chains were removed from {filename}")

    # Read PDB file
    try:
        with open(pdb, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        raise FileNotFoundError(f"PDB file not found: {pdb}")
        
    # Get available chains
    available_chains = get_available_chains(lines)
    
    # If no chains found or only one chain present
    if not available_chains:
        print(f"No chains found in {filename}")
        return
    elif len(available_chains) == 1:
        print(f"Only one chain present in {filename}: {available_chains[0]}")
        return
        
    # Handle specified keep_chain
    if keep_chain is not None:
        # Handle both string and list input for keep_chain
        if isinstance(keep_chain, str):
            keep_chain = keep_chain.upper()
            if keep_chain not in available_chains:
                raise ValueError(
                    f"Specified chain '{keep_chain}' not found in {filename}. "
                    f"Available chains: {', '.join(available_chains)}"
                )
        elif isinstance(keep_chain, (list, tuple)):
            keep_chain = [c.upper() for c in keep_chain]
            invalid_chains = [c for c in keep_chain if c not in available_chains]
            if invalid_chains:
                raise ValueError(
                    f"Specified chain(s) {', '.join(invalid_chains)} not found in {filename}. "
                    f"Available chains: {', '.join(available_chains)}"
                )
        else:
            raise TypeError("keep_chain must be a string, list, or tuple")
        write_filtered_pdb(lines, keep_chain)
        return
        
    # Interactive chain selection
    print(f"Multiple chains found in {filename}:")
    for i, chain in enumerate(available_chains):
        print(f"  {i}: {chain}")
    print("It is recommended that you remove extraneous chains.")
    
    while True:
        keep = input('Which chain(s) would you like to keep? (e.g. "A, B, C" or "all"): ').strip()
        
        if keep.lower() == "all":
            print(f"No modifications will be made to {filename}")
            return
            
        # Parse input into list of chains
        keep_chains = [c.strip().upper() for c in keep.replace(" ", "").split(",")]
        
        # Validate chains
        if all(chain in available_chains for chain in keep_chains):
            break
        print("Invalid chain ID(s). Please try again.")
    
    write_filtered_pdb(lines, keep_chains)
    
    # Ask about using same chain ID for all structures
    while True:
        same_chain = input("Would you like to use the same chain ID for all structures? (y/n): ").lower()
        if same_chain in ["y", "n"]:
            break
        print("Invalid input. Please enter 'y' or 'n'.")
    
    return keep_chains if same_chain == "y" else None

def remove_duplicates(pdb):
    """
    Removes duplicate atoms from a PDB file;
    keeping the first instance of each atom and remove the rest.
    """
    new_lines = []
    atoms = []
    with open(pdb) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_name = line[17:20]
                res_num = line[22:27]
                atom_name = line[12:16]
                if [res_name, res_num, atom_name] not in atoms:
                    atoms.append([res_name, res_num, atom_name])
                    new_lines.append(line)
            else:
                new_lines.append(line)
    with open(pdb, "w") as f:
        for line in new_lines:
            f.write(line)


def trim_pdb(ref_modified, pwd):
    """
    Trims the PDB files based on their sequence alignments with the ref
    strcuture and retains a buffer of 10 residues at each terminus.
    """
    with open(os.path.join(pwd, ref_modified)) as f:
        ref_lines = f.readlines()
    # Get the first and last residue numbers from the reference PDB file
    for line in ref_lines:
        if line.startswith("ATOM"):
            first_res = line[22:27]
            break
    for line in reversed(ref_lines):
        if line.startswith("ATOM"):
            last_res = line[22:27]
            break
    # Buffer the first and last residue numbers by 10 residues
    first_res = int(first_res) - 10
    last_res = int(last_res) + 10
    pdb_files = [i for i in os.listdir(pwd) if i.endswith("modified.pdb") and i != ref_modified]
    # Trim the PDB files
    for pdb in pdb_files:
        with open(pdb) as f:
            lines = f.readlines()
        new_lines = []
        for line in lines:
            if line.startswith("ATOM"):
                res_num = int(line[22:27])
                if res_num >= first_res and res_num <= last_res:
                    new_lines.append(line)
            else:
                new_lines.append(line)
        with open(pdb, "w") as f:
            for line in new_lines:
                f.write(line)

def multiple_prot_align():
    """
    Main funciton
    """
    check_tools(["muscle"])
    # Get the path to the directory containing the PDB files
    pwd = input("Enter the path to the directory containing the PDB files: ")

    while not os.path.exists(pwd):
        print(
            "Error: This directory does not exist. "
            "Make sure that the directory path is correct.\n"
        )
        pwd = input("Enter your working directory: ")

    os.chdir(pwd)

    # Check if any file ends with .ent
    if any(file.endswith(".ent") for file in os.listdir(pwd)):
        print('\nThe ".ent" file extensions will be updated to ".pdb".')
        for pdb_file in os.listdir(pwd):
            if pdb_file.endswith(".ent"):
                new_name = pdb_file.replace(".ent", ".pdb")
                os.rename(pdb_file, new_name)
    for pdb_file in os.listdir(pwd):
        if pdb_file.startswith("pdb"):
            new_name = pdb_file[3:]
            os.rename(pdb_file, new_name)

    # Save copies of the PDB files
    print(
        
        "Saving copies of the original PDB files to 'unaligned_structures' folder...\n"
    )
    unaligned_structures = os.path.join(pwd, "unaligned_structures")
    if not os.path.exists(unaligned_structures):
        os.mkdir(unaligned_structures)
    for pdb_file in os.listdir():
        if pdb_file.endswith(".pdb"):
            shutil.copy(pdb_file, unaligned_structures)

    # Get the name of the reference structure
    ref = input("Enter the name of the reference structure to perform an alignment: ")

    # If the user does not include the file extension, add it
    if not ref.endswith(".pdb"):
        ref += ".pdb"
    # Check that the pdb file exists
    while not os.path.isfile(ref):
        print(
            "The pdb file you entered does not appear to exist. "
            "The entry is case sensitive."
        )
        ref = input(
            "Enter the name of the reference structure to perform an alignment: "
        )
        if not ref.endswith(".pdb"):
            ref += ".pdb"

    # Run the remove_chains function
    print("Checking for extraneous chains in the PDB files...\n")
    keep_chain = None
    for pdb_file in os.listdir(pwd):
        if pdb_file.endswith(".pdb"):
            if keep_chain is None:
                keep_chain = remove_chains(os.path.join(pwd, pdb_file))
            if keep_chain:
                remove_chains(os.path.join(pwd, pdb_file), keep_chain=keep_chain)

    # Run the remove_duplicates function
    rmv_dups = input("Would you like to remove duplicate atoms from the PDB files (y/n)?: ").lower()
    while rmv_dups not in ["y", "n"]:
        print("Invalid input. Please enter 'y' or 'n'.")
        rmv_dups = input("Would you like to remove duplicate atoms from the PDB files (y/n)?: ").lower()

    if rmv_dups.lower() == "y":
        print("Removing any duplicate atoms from the PDB files...\n")
        for pdb_file in os.listdir(pwd):
            if pdb_file.endswith(".pdb"):
                remove_duplicates(os.path.join(pwd, pdb_file))
    else:
        print("Duplicate atoms will not be removed.")

    # Run the pdb_dir_to_fasta_dict function to create a dictionary of FASTA sequences
    print("Creating a dictionary of FASTA sequences...\n")
    fasta_dict = pdb_dir_to_fasta_dict(pwd)

    # Run the align_sequences function to create a dictionary of aligned sequences
    print("Aligning the sequences...\n")
    aligned_dict = align_sequences(fasta_dict, "aligned.fasta")

    # Run the get_aligned_residues function to create a dictionary of
    # new residue positions
    print("Determining the new residue positions based on the alignment...\n")
    new_positions = get_aligned_residues(aligned_dict)

    # Run the modify_pdb_sequence function
    print("Modifying the PDB files with new residue positions...\n")
    for pdb_file in os.listdir(pwd):
        if pdb_file.endswith(".pdb"):
            if pdb_file[:-4] in new_positions:
                modify_pdb_sequence(
                    pdb_file,
                    new_positions[pdb_file[:-4]],
                    f"{pdb_file[:-4]}_modified.pdb",
                )
        else:
            continue
        # Check that the modified PDB file was created
        if os.path.exists(f"{pdb_file[:-4]}_modified.pdb"):
            print(f"The modified PDB file for {pdb_file} was created successfully.")
        else:
            print(f"The modified PDB file for {pdb_file} was NOT created successfully.")

    # Run the get_og_pos function
    print(
        "Determining the original residue positions " "for each modified PDB file...\n"
    )
    pos_dfs = {}
    for pdb_file in os.listdir(pwd):
        if pdb_file.endswith(".pdb"):
            if pdb_file[:-4] in new_positions:
                df = get_og_pos(
                    pdb_file, f"{os.path.splitext(pdb_file)[0]}_modified.pdb"
                )
                pos_dfs[pdb_file[:-4]] = df
        else:
            continue

    # Merge the dataframes
    merged_df = pos_dfs[list(pos_dfs.keys())[0]]
    for key in pos_dfs:
        if key == list(pos_dfs.keys())[0]:
            continue
        merged_df = pd.merge(merged_df, pos_dfs[key], how="outer", on=["new_pos"])
    merged_df["new_pos"] = merged_df["new_pos"].astype(int)
    merged_df = merged_df.sort_values(by=["new_pos"])
    # Write the merged dataframe to a CSV file
    print("Writing the new and original residue positions to a CSV file...\n")
    merged_df.to_csv("residue_positions.csv", index=False)

    # Run the trim_pdb function
    print("Trimming the PDB files with respect to the reference...\n")
    trim_pdb(f"{ref[:-4]}_modified.pdb", pwd)

    # Run the align_structures function
    print("Aligning structures to the reference...")
    align_structures(pwd, f"{ref[:-4]}_modified.pdb")

    # Move all non-structural files to a new folder called "results"
    print('\nMoving all non-structural files to a new folder called "results"...')
    if not os.path.exists("results"):
        os.mkdir("results")
    else:
        print(
            'The "results" folder already exists. Files will be moved to this folder.'
        )
    for file in os.listdir(pwd):
        if not file.endswith(".pdb") and not os.path.isdir(file):
            try:
                shutil.move(file, "results")
            except:
                print(
                    f'Unable to move {file} to the "results" folder. '
                    "It may already exist in the folder."
                )

    # Remove all of the PDB files from the working directory that
    # do not end with "_modified.pdb"
    dir_files = [i for i in os.listdir() if i.endswith(".pdb") and not i.endswith("_modified.pdb")]
    for file in dir_files:
        os.remove(file)
    print("DONE!")


if __name__ == "__main__":
    multiple_prot_align()
