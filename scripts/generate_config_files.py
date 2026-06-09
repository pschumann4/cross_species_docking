import os
import re
import numpy as np
from collections import defaultdict
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import euclidean3d, centroid


def calculate_gridbox_from_pdb(pdb_file, ligand_resname):
    """
    Calculate gridbox parameters directly from a reference PDB file containing a ligand.
    
    METHODOLOGY:
    ------------
    1. Extract all HETATM coordinates for the specified ligand
    2. Calculate the ligand centroid (geometric center of all ligand atoms)
    3. Find the maximum distance from centroid to any ligand atom (ligand radius)
    4. Set grid box size as 4× this maximum radius in all dimensions (cubic)
    5. Center the grid box at the ligand centroid coordinates

    SCIENTIFIC RATIONALE:
    ---------------------
    The 4× multiplier ensures the grid box fully encompasses the ligand with adequate
    padding for rotational and translational sampling during docking. This approach:
    - Captures the native binding pose completely
    - Allows conformational flexibility during docking
    - Provides consistent grid sizing across different ligand geometries
    - Produces a cubic box guaranteed to be large enough in every direction
    
    Parameters:
    -----------
    pdb_file : str
        Path to PDB file containing the protein-ligand complex
    ligand_resname : str
        Three-letter residue name of the ligand (e.g., 'ATP', 'HEM')
    
    Returns:
    --------
    dict : Gridbox parameters with keys:
        'size_x', 'size_y', 'size_z' : Grid dimensions (Angstroms, all equal)
        'center_x', 'center_y', 'center_z' : Grid center coordinates
    None : If ligand not found or error occurs
    """
    ligand_coords = []
    
    # Extract ligand coordinates from PDB file
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("HETATM"):
                    # PDB format: residue name at columns 17-20 (0-indexed: 17:20)
                    if line[17:20].strip() == ligand_resname:
                        # PDB format: x at 30-38, y at 38-46, z at 46-54
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append([x, y, z])
    except Exception as e:
        print(f"Error reading PDB file {pdb_file}: {e}")
        return None
    
    # Check if ligand was found
    if not ligand_coords:
        print(f"Warning: Ligand '{ligand_resname}' not found in {pdb_file}")
        return None
    
    print(f"Found {len(ligand_coords)} atoms for ligand '{ligand_resname}'")

    coords_arr = np.array(ligand_coords)

    # Calculate ligand centroid (geometric centre)
    ligand_centroid = centroid(ligand_coords)
    print(f"Ligand centroid: [{ligand_centroid[0]:.3f}, {ligand_centroid[1]:.3f}, {ligand_centroid[2]:.3f}]")

    # Calculate grid box dimensions
    max_radius = max(euclidean3d(ligand_centroid, coord) for coord in ligand_coords)
    size = round(4 * max_radius, 3)
    size_x = size_y = size_z = size

    print(f"Max ligand radius: {max_radius:.3f} Å")
    print(f"Grid box size:     {size_x:.3f} × {size_y:.3f} × {size_z:.3f} Å  (4 × max_radius, cubic)")

    # Create gridbox parameters dictionary
    gridbox_params = {
        'size_x': str(size_x),
        'size_y': str(size_y),
        'size_z': str(size_z),
        'center_x': str(round(ligand_centroid[0], 3)),
        'center_y': str(round(ligand_centroid[1], 3)),
        'center_z': str(round(ligand_centroid[2], 3))
    }
    
    return gridbox_params


def read_gridbox_file(filepath):
    """
    Read gridbox parameters from a text file.
    Returns a dictionary containing the gridbox parameters.
    
    Expected format:
    PDB file: <pdb_name>
    size_x: <value>
    size_y: <value>
    size_z: <value>
    center_x: <value>
    center_y: <value>
    center_z: <value>
    """
    gridbox_params = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if ':' in line:
                    key, value = line.strip().split(':', 1)  # Split only on first colon
                    key = key.strip().lower()
                    value = value.strip()
                    if key != 'pdb file':
                        # Remove any possible unit or additional text and convert to string
                        value = str(float(value.split()[0]))
                        gridbox_params[key] = value
        
        # Verify all required parameters are present
        required_params = ['size_x', 'size_y', 'size_z', 'center_x', 'center_y', 'center_z']
        for param in required_params:
            if param not in gridbox_params:
                raise ValueError(f"Missing required parameter: {param}")
        
        return gridbox_params
    except Exception as e:
        print(f"Error reading gridbox file: {str(e)}")
        return None


def parse_ensemble_filename(filename):
    """
    Parse ensemble filenames to extract base name and ensemble number.
    
    Examples:
        'protein_ensemble_1.pdbqt' -> ('protein', 1)
        'protein_ensemble_1_rigid.pdbqt' -> ('protein', 1)
        'protein_ensemble_1_flex.pdbqt' -> ('protein', 1)
        'protein.pdbqt' -> ('protein', None)
    
    Returns:
        tuple: (base_name, ensemble_number or None)
    """
    # Remove .pdbqt extension
    name = filename.replace('.pdbqt', '')
    
    # Remove _rigid or _flex suffix if present
    name = name.replace('_rigid', '').replace('_flex', '')
    
    # Check for ensemble pattern: _ensemble_N
    ensemble_pattern = r'(.+)_ensemble_(\d+)$'
    match = re.match(ensemble_pattern, name)
    
    if match:
        base_name = match.group(1)
        ensemble_num = int(match.group(2))
        return (base_name, ensemble_num)
    else:
        return (name, None)


def read_flexible_residues(residues_file):
    """
    Read flexible residues from the summary file created by mds_structure_prep.py
    
    The expected format is:
    FLEXIBLE RESIDUES SUMMARY
    ==========================================================
    
    structure_name: [resid1,resid2,resid3]
    another_structure: []
    
    ==========================================================
    Total structures analyzed: N
    Structures with flexible residues: M
    
    Parameters:
    -----------
    residues_file : str
        Path to the flex_residues.txt file
    
    Returns:
    --------
    dict: Dictionary mapping (base_name, ensemble_num) to has_flex boolean
    """
    flex_info = {}
    
    with open(residues_file, 'r') as f:
        lines = f.readlines()
    
    # Parse each line
    for line in lines:
        line = line.strip()
        
        # Skip empty lines, header, separator lines, and statistics
        if not line or line.startswith('=') or line.startswith('FLEXIBLE') or \
           line.startswith('Total') or line.startswith('Structures with'):
            continue
        
        # Parse structure_name: [resid1,resid2,resid3] format
        if ':' in line and '[' in line and ']' in line:
            try:
                # Split on first colon
                parts = line.split(':', 1)
                if len(parts) != 2:
                    continue
                
                file_identifier = parts[0].strip()
                residue_info = parts[1].strip()
                
                # Extract residue list (between brackets)
                residue_info = residue_info.split('[')[1].split(']')[0].strip()
                
                # Check if there are any flexible residues
                has_flex = bool(residue_info)
                
                # Parse the file identifier to get base name and ensemble number
                base_name, ensemble_num = parse_ensemble_filename(file_identifier)
                if ensemble_num is None:
                    ensemble_num = 0
                
                flex_info[(base_name, ensemble_num)] = has_flex
                
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line: {line}")
                print(f"  Error: {e}")
                continue
    
    return flex_info


def group_pdbqt_files(pdbqt_files):
    """
    Group PDBQT files by their base name and ensemble number.
    
    Returns:
        dict: {base_name: {ensemble_num: {'rigid': filename, 'flex': filename, 'full': filename}}}
    """
    grouped = defaultdict(lambda: defaultdict(dict))
    
    for filename in pdbqt_files:
        base_name, ensemble_num = parse_ensemble_filename(filename)
        
        # Determine file type
        if '_rigid.pdbqt' in filename:
            file_type = 'rigid'
        elif '_flex.pdbqt' in filename:
            file_type = 'flex'
        else:
            file_type = 'full'
        
        # Use ensemble_num=0 for non-ensemble files
        if ensemble_num is None:
            ensemble_num = 0
        
        grouped[base_name][ensemble_num][file_type] = filename.replace('.pdbqt', '')
    
    return grouped


def get_config_files():
    """
    Create AutoDock Vina config files for each PDBQT file, including ensemble members.
    Handles both single structures and ensemble structures with proper naming.

    Reads project_config.json to pre-fill ligand name, reference PDB path, ligand
    residue name, and docking parameters collected in earlier steps.  Any new values
    entered here are written back to the config for subsequent scripts.
    """
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from utils import (resolve_project_dir, load_config, save_config,
                       get_project_paths, resolve_reference_pdb)

    project_dir = resolve_project_dir()
    config = load_config(project_dir)
    paths = get_project_paths(project_dir)
    pdbqt_dir = paths["pdbqt_files"]

    if not os.path.isdir(pdbqt_dir):
        print(f"Error: pdbqt_files/ not found at {pdbqt_dir}")
        print("Please run prep_pdbqt.py first.")
        return
    os.chdir(pdbqt_dir)

    # ── Ligand name ──────────────────────────────────────────────────────────
    ligand_name = config.get("ligand_name")
    if ligand_name:
        if os.path.exists(f"{ligand_name}.pdbqt"):
            print(f"\nUsing ligand from config: {ligand_name}.pdbqt")
        else:
            print(f"\nWarning: {ligand_name}.pdbqt not found in pdbqt_files/. Please enter manually.")
            ligand_name = None

    if not ligand_name:
        same_dir = input("Is the ligand PDBQT file in the 'pdbqt_files' directory? (y/n): ").lower()
        if same_dir == "n":
            ligand_path = input("Enter the path to the ligand PDBQT file: ").strip('"')
            while not os.path.exists(ligand_path):
                ligand_path = input("That path does not appear to exist.\nPlease enter the path to the ligand PDBQT file: ").strip('"')
            ligand_name = os.path.basename(ligand_path).replace('.pdbqt', '')
            try:
                import shutil as _shutil
                _shutil.copy(ligand_path, pdbqt_dir)
                print(f"Copied ligand file to {pdbqt_dir}")
            except Exception as e:
                print(f"Error copying ligand file: {e}")
                return
        else:
            ligand_name = input(
                "\nEnter the name of the ligand PDBQT file that will be used for the docking simulation: "
            )
            if ligand_name.endswith(".pdbqt"):
                ligand_name = ligand_name.replace(".pdbqt", "")
        config["ligand_name"] = ligand_name
        save_config(project_dir, config)

    # ── Gridbox parameters ───────────────────────────────────────────────────
    reference_pdb = resolve_reference_pdb(config, project_dir)
    ligand_resname = config.get("ligand_resname")
    gridbox_params = None

    if reference_pdb and ligand_resname:
        print(f"\nAuto-extracting gridbox from: {os.path.basename(reference_pdb)}"
              f"  (ligand: {ligand_resname})")
        print("-"*60)
        gridbox_params = calculate_gridbox_from_pdb(reference_pdb, ligand_resname)
        if gridbox_params is None:
            print("Auto-extraction failed. Falling back to manual selection.")

    if gridbox_params is None:
        print("\n" + "="*60)
        print("GRIDBOX PARAMETER CONFIGURATION")
        print("="*60)
        print("\nChoose gridbox parameter source:")
        print("  1. Automatically extract from reference PDB file (recommended)")
        print("  2. Use existing gridbox coordinate file")

        while True:
            gridbox_choice = input("\nSelect option (1 or 2): ").strip()
            if gridbox_choice in ['1', '2']:
                break
            print("Invalid choice. Please enter 1 or 2.")

        if gridbox_choice == '1':
            print("\n" + "-"*60)
            print("AUTOMATIC GRIDBOX EXTRACTION")
            print("-"*60)

            # Scan project_dir for a ref_ PDB if not already resolved
            if reference_pdb is None:
                pdb_candidates = [f for f in os.listdir(project_dir)
                                  if f.startswith('ref_') and f.endswith('.pdb')]
                if len(pdb_candidates) == 1:
                    candidate = os.path.join(project_dir, pdb_candidates[0])
                    print(f"\nFound reference PDB: {pdb_candidates[0]}")
                    if input("Use this file? (y/n): ").lower() == 'y':
                        reference_pdb = candidate

            while reference_pdb is None:
                reference_pdb = input("\nEnter path to reference PDB file: ").strip().strip('"')
                if not os.path.exists(reference_pdb):
                    print("That path does not exist.")
                    reference_pdb = None

            if not ligand_resname:
                print("\nEnter the 3-letter residue name of the ligand in the reference PDB")
                ligand_resname = input("Ligand residue name: ").strip().upper()
                config["ligand_resname"] = ligand_resname
                save_config(project_dir, config)

            print("\nExtracting gridbox parameters...")
            print("-"*60)
            gridbox_params = calculate_gridbox_from_pdb(reference_pdb, ligand_resname)

            if gridbox_params is None:
                print("\nFailed to extract gridbox parameters.")
                print("Falling back to existing gridbox file option...")
                gridbox_choice = '2'

        if gridbox_choice == '2' or gridbox_params is None:
            gridbox_file = None
            details_dir = paths["results"]
            if os.path.isdir(details_dir):
                gridbox_files = [f for f in os.listdir(details_dir) if f.endswith('_gridbox_coords.txt')]
                if len(gridbox_files) == 1:
                    gridbox_file = os.path.join(details_dir, gridbox_files[0])
                    print(f"\nFound gridbox file in results/: {gridbox_files[0]}")
                elif len(gridbox_files) > 1:
                    print(f"\nFound multiple gridbox files in results/:")
                    for i, f in enumerate(gridbox_files, 1):
                        print(f"  {i}. {f}")
                    while True:
                        choice = input(f"Select file (1-{len(gridbox_files)}): ")
                        if choice.isdigit() and 1 <= int(choice) <= len(gridbox_files):
                            gridbox_file = os.path.join(details_dir, gridbox_files[int(choice) - 1])
                            break
                        print("Invalid selection.")

            while True:
                if gridbox_file is None:
                    gridbox_file = input("Enter the path to the gridbox parameters file: ").strip().strip('"')
                else:
                    print(f"Reading gridbox parameters from: {gridbox_file}")
                if not os.path.exists(gridbox_file):
                    print("That path does not appear to exist.")
                    gridbox_file = None
                    continue
                gridbox_params = read_gridbox_file(gridbox_file)
                if gridbox_params:
                    break
                gridbox_file = None

    if gridbox_params is None:
        print("\nError: Failed to obtain gridbox parameters. Exiting.")
        return

    print("\n" + "="*60)
    print("GRIDBOX PARAMETERS")
    print("="*60)
    print(f"Grid box size:   {gridbox_params['size_x']} × {gridbox_params['size_y']} × {gridbox_params['size_z']} Å")
    print(f"Grid box center: ({gridbox_params['center_x']}, {gridbox_params['center_y']}, {gridbox_params['center_z']})")

    # Get all PDBQT files in the directory
    pdbqt_files = [f for f in os.listdir(pdbqt_dir) if f.endswith(".pdbqt") and not f.startswith(ligand_name)]

    if not pdbqt_files:
        print("No PDBQT files found in the specified directory.")
        return

    # Group PDBQT files by base name and ensemble number
    grouped_files = group_pdbqt_files(pdbqt_files)

    # Display summary of files
    print("\n" + "="*60)
    print("PDBQT FILES DETECTED")
    print("="*60)

    single_count = 0
    ensemble_count = 0

    for base_name in sorted(grouped_files.keys()):
        ensemble_nums = list(grouped_files[base_name].keys())
        if 0 in ensemble_nums and len(ensemble_nums) == 1:
            single_count += 1
            print(f"\nSingle: {base_name}")
        else:
            actual_ensemble_nums = [n for n in ensemble_nums if n != 0]
            if actual_ensemble_nums:
                ensemble_count += len(actual_ensemble_nums)
                print(f"\nEnsemble: {base_name}")
                print(f"  Members: {len(actual_ensemble_nums)} structures (ensemble_{min(actual_ensemble_nums)} to ensemble_{max(actual_ensemble_nums)})")

    print(f"\nTotal: {single_count} single structures, {ensemble_count} ensemble members")
    print("="*60)

    # ── Docking parameters (pre-filled from config, press Enter to keep) ─────
    saved_docking = config.get("docking", {})
    if saved_docking:
        print("\nDocking parameters from config (press Enter to keep each value):")

    while True:
        scoring = input(f"Scoring function (ad4 / vina / vinardo) [{saved_docking.get('scoring', 'vina')}]: ").strip()
        if not scoring:
            scoring = saved_docking.get("scoring", "vina")

        num_modes = input(f"Number of modes (1–20) [{saved_docking.get('num_modes', '5')}]: ").strip()
        if not num_modes:
            num_modes = saved_docking.get("num_modes", "5")

        energy_range = input(f"Energy range kcal/mol (1–10) [{saved_docking.get('energy_range', '3')}]: ").strip()
        if not energy_range:
            energy_range = saved_docking.get("energy_range", "3")

        exhaustiveness = input(f"Exhaustiveness (8–32) [{saved_docking.get('exhaustiveness', '16')}]: ").strip()
        if not exhaustiveness:
            exhaustiveness = saved_docking.get("exhaustiveness", "16")

        print("\n" + "="*60)
        print("CONFIGURATION PARAMETERS")
        print("="*60)
        print(f"Ligand: {ligand_name}.pdbqt")
        print(f"Scoring function: {scoring}")
        print(f"Number of modes: {num_modes}")
        print(f"Energy range: {energy_range}")
        print(f"Exhaustiveness: {exhaustiveness}")
        print(f"Grid box size: {gridbox_params['size_x']} × {gridbox_params['size_y']} × {gridbox_params['size_z']} Å")
        print(f"Grid box center: ({gridbox_params['center_x']}, {gridbox_params['center_y']}, {gridbox_params['center_z']})")
        print("="*60)

        confirm = input("\nIs this information correct? (y/n): ").lower()
        while confirm not in ["y", "n"]:
            print("Invalid input.")
            confirm = input("Is this information correct? (y/n): ").lower()

        if confirm == "y":
            config["docking"] = {
                "scoring": scoring,
                "num_modes": num_modes,
                "energy_range": energy_range,
                "exhaustiveness": exhaustiveness,
            }
            save_config(project_dir, config)
            break

    # ── Flexible residues (auto-detected from results/flex_residues.txt) ──────
    flex_info = {}
    details_residues_file = os.path.join(paths["results"], "flex_residues.txt")

    if os.path.exists(details_residues_file):
        print(f"\nFound flex_residues.txt in results/ — loading flexible residue data...")
        try:
            flex_info = read_flexible_residues(details_residues_file)
            structures_with_flex = sum(1 for has_flex in flex_info.values() if has_flex)
            print(f"  Loaded data for {len(flex_info)} structure(s): "
                  f"{structures_with_flex} flexible, "
                  f"{len(flex_info) - structures_with_flex} rigid")
        except Exception as e:
            print(f"  Warning: could not read flex_residues.txt ({e}) — treating all receptors as rigid.")
            flex_info = {}
    else:
        print("\nNo flex_residues.txt found in results/ — all receptors will be treated as rigid.")

    # Categorize files
    flexible_configs = []
    rigid_configs = []

    for base_name in grouped_files:
        for ensemble_num in grouped_files[base_name]:
            file_info = grouped_files[base_name][ensemble_num]

            if base_name == ligand_name:
                continue

            has_flex = flex_info.get((base_name, ensemble_num), False)
            if not has_flex and ensemble_num != 0:
                has_flex = flex_info.get((base_name, 0), False)

            if has_flex:
                if 'rigid' in file_info and 'flex' in file_info:
                    flexible_configs.append((base_name, ensemble_num, file_info))
                else:
                    label = base_name if ensemble_num == 0 else f"{base_name}_ensemble_{ensemble_num}"
                    print(f"Warning: {label} marked as flexible but missing rigid/flex files — skipping")
            else:
                if 'full' in file_info:
                    rigid_configs.append((base_name, ensemble_num, file_info))
                elif 'rigid' in file_info:
                    rigid_configs.append((base_name, ensemble_num, file_info))

    # Display categorization
    print("\n" + "="*60)
    print("FILE CATEGORIZATION")
    print("="*60)
    print("\nFLEXIBLE RECEPTORS:")
    if flexible_configs:
        for base_name, ens_num, _ in flexible_configs:
            print(f"  {base_name}" if ens_num == 0 else f"  {base_name}_ensemble_{ens_num}")
    else:
        print("  (none)")

    print("\nRIGID RECEPTORS:")
    if rigid_configs:
        for base_name, ens_num, _ in rigid_configs:
            print(f"  {base_name}" if ens_num == 0 else f"  {base_name}_ensemble_{ens_num}")
    else:
        print("  (none)")
    print("="*60)
    
    # Create config files
    config_count = 0
    
    # Flexible receptors
    for base_name, ens_num, file_info in flexible_configs:
        # Generate config filename
        if ens_num == 0:
            config_name = f"{base_name}_conf.txt"
            output_name = f"{base_name}_bound_{ligand_name}.pdbqt"
        else:
            config_name = f"{base_name}_ensemble_{ens_num}_conf.txt"
            output_name = f"{base_name}_ensemble_{ens_num}_bound_{ligand_name}.pdbqt"
        
        with open(config_name, "w") as f:
            f.write(f"flex = {file_info['flex']}.pdbqt\n")
            f.write(f"receptor = {file_info['rigid']}.pdbqt\n")
            f.write(f"ligand = {ligand_name}.pdbqt\n")
            f.write(f"scoring = {scoring}\n\n")
            f.write(f"size_x = {gridbox_params['size_x']}\n")
            f.write(f"size_y = {gridbox_params['size_y']}\n")
            f.write(f"size_z = {gridbox_params['size_z']}\n\n")
            f.write(f"center_x = {gridbox_params['center_x']}\n")
            f.write(f"center_y = {gridbox_params['center_y']}\n")
            f.write(f"center_z = {gridbox_params['center_z']}\n\n")
            f.write(f"spacing = 1\n\n")
            f.write(f"num_modes = {num_modes}\n")
            f.write(f"energy_range = {energy_range}\n")
            f.write(f"exhaustiveness = {exhaustiveness}\n\n")
            f.write(f"out = {output_name}")
        
        config_count += 1
    
    # Rigid receptors
    for base_name, ens_num, file_info in rigid_configs:
        # Generate config filename
        if ens_num == 0:
            config_name = f"{base_name}_conf.txt"
            output_name = f"{base_name}_bound_{ligand_name}.pdbqt"
        else:
            config_name = f"{base_name}_ensemble_{ens_num}_conf.txt"
            output_name = f"{base_name}_ensemble_{ens_num}_bound_{ligand_name}.pdbqt"
        
        # Use 'full' if available, otherwise 'rigid'
        receptor_file = file_info.get('full', file_info.get('rigid'))
        
        with open(config_name, "w") as f:
            f.write(f"receptor = {receptor_file}.pdbqt\n")
            f.write(f"ligand = {ligand_name}.pdbqt\n")
            f.write(f"scoring = {scoring}\n\n")
            f.write(f"size_x = {gridbox_params['size_x']}\n")
            f.write(f"size_y = {gridbox_params['size_y']}\n")
            f.write(f"size_z = {gridbox_params['size_z']}\n\n")
            f.write(f"center_x = {gridbox_params['center_x']}\n")
            f.write(f"center_y = {gridbox_params['center_y']}\n")
            f.write(f"center_z = {gridbox_params['center_z']}\n\n")
            f.write(f"spacing = 1\n\n")
            f.write(f"num_modes = {num_modes}\n")
            f.write(f"energy_range = {energy_range}\n")
            f.write(f"exhaustiveness = {exhaustiveness}\n\n")
            f.write(f"out = {output_name}")
        
        config_count += 1
    
    print(f"\n{'='*60}")
    print(f"SUCCESS: Created {config_count} configuration files")
    print(f"{'='*60}")
    print(f"Files saved to: {os.path.abspath(pdbqt_dir)}")


if __name__ == "__main__":
    get_config_files()