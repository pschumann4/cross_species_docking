import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import parse_hetatm_coords, rmsd_hungarian


def ligand_rmsd():
    """
    Calculate the RMSD between a reference ligand and the docked ligand in each
    model PDB file, then write a summary text file and histogram.

    The best-RMSD pose selection that was previously done here has been moved
    upstream into run_vina_batch.py (process_vina_output_top_mode_only).  Each
    model PDB on disk already contains the pose most similar to the reference
    binding mode, so this script is now a clean reporting tool: it reads the
    saved model, computes one RMSD per species, and records the result.

    Input:
    - Directory containing *_model.pdb files and a reference PDB.
    - The ligand residue name (three-letter code) as it appears in the PDB files.

    Output:
    - ligand_rmsd.txt  — one RMSD value per model, written to the nearest
                         'results' directory (or the working directory if none found).
    - ligand_rmsd.png  — histogram of RMSD values.
    """
    from utils import (resolve_project_dir, load_config, save_config,
                       get_project_paths, resolve_reference_pdb)

    project_dir = resolve_project_dir()
    config = load_config(project_dir)
    paths = get_project_paths(project_dir)

    # ── Models directory ──────────────────────────────────────────────────────
    pdb_dir = paths["models"]
    if not os.path.isdir(pdb_dir):
        print(f"\nWarning: models/ not found at {pdb_dir}")
        pdb_dir = input("Enter the directory containing the PDB models: ").strip().strip('"')
        while not os.path.exists(pdb_dir):
            print("The working directory does not exist.")
            pdb_dir = input("Enter the directory: ").strip().strip('"')
    else:
        print(f"\nUsing models directory: {pdb_dir}")

    os.chdir(pdb_dir)

    # Output always goes to the project results/ directory
    output_dir = paths["results"]
    os.makedirs(output_dir, exist_ok=True)

    # ── Locate reference PDB ──────────────────────────────────────────────────
    pdb_files = [i for i in os.listdir(pdb_dir) if i.endswith(".pdb")]
    ref_pdb = ""

    # Try the config first (the reference is copied into models/ during run_vina_batch)
    ref_basename = config.get("reference_pdb")
    if ref_basename and os.path.exists(os.path.join(pdb_dir, ref_basename)):
        ref_pdb = ref_basename
        print(f"Using reference PDB from config: {ref_pdb}")
    else:
        # Fall back to interactive scan for ref_ files
        ref = "n"
        for file in pdb_files:
            if file.startswith("ref_"):
                ref = input(f"Is {file} the reference PDB file? (y/n): ").lower()
                while ref not in ["y", "n"]:
                    ref = input("Please enter y or n: ")
                if ref == "y":
                    ref_pdb = file
                    break
        if ref == "n":
            ref_pdb = input("Enter the name of the reference PDB file (if none, type 'random'): ")
            if not ref_pdb.endswith(".pdb"):
                ref_pdb += ".pdb"

    while not os.path.exists(ref_pdb):
        print("The reference file does not exist.")
        ref_pdb = input("Enter the name of the reference PDB file: ")
        if not ref_pdb.endswith(".pdb"):
            ref_pdb += ".pdb"

    # ── Ligand residue name ───────────────────────────────────────────────────
    ligand = config.get("ligand_resname")
    if ligand:
        print(f"Using ligand ID from config: {ligand}")
    else:
        ligand = input("Enter the ligand ID as it is found within the PDB files: ")

    # Extract reference ligand heavy-atom coords
    with open(ref_pdb, "r") as f:
        ref_lines = f.readlines()

    coords1, elements1 = parse_hetatm_coords(ref_lines, ligand)

    while len(coords1) == 0:
        print("ERROR: The ligand entered was not found in the PDB file.")
        ligands = {l[17:20].strip() for l in ref_lines if l.startswith("HETATM")}
        print(f"The following ligands were found in the PDB file: {ligands}")
        ligand = input("Please select one of the ligands from the list: ")
        coords1, elements1 = parse_hetatm_coords(ref_lines, ligand)

    # Calculate RMSD for every model PDB
    rmsd_dict = {}
    model_files = [
        i for i in os.listdir(pdb_dir)
        if i.endswith(".pdb") and "model" in i and i != ref_pdb
    ]

    for file in model_files:
        base = os.path.splitext(os.path.basename(file))[0]

        with open(file, "r") as f:
            model_lines = f.readlines()

        coords2, elements2 = parse_hetatm_coords(model_lines, ligand)

        if len(coords2) == 0:
            print(f"The ligand was not found in {file}.")
            continue

        if len(elements1) != len(elements2):
            print(f"ERROR: Different number of atoms in {ref_pdb} "
                  f"({len(elements1)}) and {file} ({len(elements2)})")
            continue

        if sorted(elements1) != sorted(elements2):
            print(f"ERROR: Different element types in {ref_pdb} and {file}")
            print(f"  {ref_pdb}: {sorted(elements1)}")
            print(f"  {file}: {sorted(elements2)}")
            continue

        rmsd_dict[base] = round(rmsd_hungarian(coords1, coords2), 3)

    # Write results
    rmsd_txt_path = os.path.join(output_dir, "ligand_rmsd.txt")
    with open(rmsd_txt_path, "w") as f:
        f.write("Ligand RMSD values\n\n")
        for name, rmsd in rmsd_dict.items():
            f.write(f"{name}: {rmsd}\n")

    if os.path.exists(rmsd_txt_path):
        print(f"The ligand RMSD values were written to {rmsd_txt_path}.")

    # Histogram
    rmsd_values = list(rmsd_dict.values())
    plt.hist(rmsd_values, bins=30)
    plt.xlabel("RMSD (Å)", fontsize=18)
    plt.ylabel("Frequency", fontsize=18)
    for patch in plt.gca().patches:
        patch.set_facecolor("#A4D4F7")
        patch.set_edgecolor("black")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "ligand_rmsd.png"), dpi=300)
    plt.show()


if __name__ == "__main__":
    ligand_rmsd()
