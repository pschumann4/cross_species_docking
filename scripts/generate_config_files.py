import os
import re
import numpy as np
from collections import defaultdict


def euclidean3d(v1, v2):
    """
    Calculate euclidean distance between two 3D points.
    Integrated from get_pocket.py
    """
    if not (len(v1) == 3 and len(v2) == 3):
        return None
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def centroid(coords):
    """
    Calculate the centroid from a 3D point cloud.
    Integrated from get_pocket.py
    
    Parameters:
    -----------
    coords : list of lists
        Array of [x, y, z] coordinate arrays
    
    Returns:
    --------
    list : [x, y, z] centroid coordinates
    """
    return list(
        map(
            np.mean,
            (([c[0] for c in coords]), ([c[1] for c in coords]), ([c[2] for c in coords])),
        )
    )


def calculate_gridbox_from_pdb(pdb_file, ligand_resname):
    """
    Calculate gridbox parameters directly from a reference PDB file containing a ligand.
    
    METHODOLOGY (from get_pocket.py):
    ---------------------------------
    1. Extract all HETATM coordinates for the specified ligand
    2. Calculate the ligand centroid (geometric center of all ligand atoms)
    3. Find the maximum distance from centroid to any ligand atom (ligand radius)
    4. Set grid box size as 4× this maximum radius in all dimensions
    5. Center the grid box at the ligand centroid coordinates
    
    SCIENTIFIC RATIONALE:
    --------------------
    The 4× multiplier ensures the grid box fully encompasses the ligand with adequate
    padding for rotational and translational sampling during docking. This approach:
    - Captures the native binding pose completely
    - Allows conformational flexibility during docking
    - Provides consistent grid sizing across different ligand geometries
    - Follows the PLIP-based binding site definition methodology
    
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
    
    # Calculate ligand centroid
    ligand_centroid = centroid(ligand_coords)
    print(f"Ligand centroid: [{ligand_centroid[0]:.3f}, {ligand_centroid[1]:.3f}, {ligand_centroid[2]:.3f}]")
    
    # Calculate maximum distance from centroid to any ligand atom
    max_radius = 0
    for coords in ligand_coords:
        dist = euclidean3d(ligand_centroid, coords)
        if dist > max_radius:
            max_radius = dist
    
    print(f"Maximum ligand radius: {max_radius:.3f} Å")
    
    # Calculate grid box size: 4× maximum radius ensures full coverage
    # This multiplier provides adequate space for ligand rotation and translation
    grid_size = round(max_radius * 4, 3)
    
    print(f"Calculated grid box size: {grid_size:.3f} Å (cubic)")
    
    # Create gridbox parameters dictionary
    gridbox_params = {
        'size_x': str(grid_size),
        'size_y': str(grid_size),
        'size_z': str(grid_size),
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
    
    NEW FEATURE: Can automatically extract gridbox parameters from a reference PDB file
    containing a protein-ligand complex, eliminating the need for pre-generated gridbox
    coordinate files.
    """
    # Ask the user for the directory where the PDBQT files are located
    pwd = input("Enter the path to the 'pdbqt_files' directory: ")

    # Check if the path exists
    while not os.path.exists(pwd):
        pwd = input(
            "That path does not appear to exist.\nPlease enter the path to "
            "the directory containing the PDBQT files: "
        )
    os.chdir(pwd)

    # Ask the user for the name of the ligand PDBQT file that will be used for the docking simulation
    ligand_name = input(
        "\nEnter the name of the ligand PDBQT file that will be used for the docking simulation: "
    )
    # If the ligand name ends with .pdbqt, remove it
    if ligand_name.endswith(".pdbqt"):
        ligand_name = ligand_name.replace(".pdbqt", "")

    # Gridbox parameter acquisition
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
    
    gridbox_params = None
    
    if gridbox_choice == '1':
        # AUTOMATIC EXTRACTION FROM REFERENCE PDB
        print("\n" + "-"*60)
        print("AUTOMATIC GRIDBOX EXTRACTION")
        print("-"*60)
        print("\nThis will extract binding pocket coordinates from a reference")
        print("PDB file containing a protein-ligand complex.")
        
        # Try to find reference PDB automatically
        parent_dir = os.path.dirname(pwd)
        reference_pdb = None
        
        # Look for PDB files in parent directory
        if os.path.exists(parent_dir):
            pdb_files = [f for f in os.listdir(parent_dir) 
                        if f.startswith('ref_') and f.endswith('.pdb')]
            
            if len(pdb_files) == 1:
                reference_pdb = os.path.join(parent_dir, pdb_files[0])
                print(f"\nFound reference PDB in parent directory: {pdb_files[0]}")
                use_found = input("Use this file? (y/n): ").lower()
                if use_found != 'y':
                    reference_pdb = None

        # Prompt for path if not found automatically
        while reference_pdb is None:
            reference_pdb = input("\nEnter path to reference PDB file: ").strip()
            if reference_pdb.startswith('"') and reference_pdb.endswith('"'):
                reference_pdb = reference_pdb[1:-1]
            
            if not os.path.exists(reference_pdb):
                print("That path does not exist.")
                reference_pdb = None
        
        # Get ligand residue name
        print("\nEnter the 3-letter residue name of the ligand in the reference PDB")
        ligand_resname = input("Ligand residue name: ").strip().upper()
        
        # Extract gridbox parameters
        print("\nExtracting gridbox parameters...")
        print("-"*60)
        gridbox_params = calculate_gridbox_from_pdb(reference_pdb, ligand_resname)
        
        if gridbox_params is None:
            print("\nFailed to extract gridbox parameters.")
            print("Falling back to existing gridbox file option...")
            gridbox_choice = '2'
    
    if gridbox_choice == '2' or gridbox_params is None:
        # EXISTING GRIDBOX FILE
        parent_dir = os.path.dirname(pwd)
        gridbox_file = None
        
        # Search in ../details/ for files ending with _gridbox_coords.txt
        details_dir = os.path.join(parent_dir, "details")
        if os.path.exists(details_dir):
            gridbox_files = [f for f in os.listdir(details_dir) if f.endswith('_gridbox_coords.txt')]
            if len(gridbox_files) == 1:
                gridbox_file = os.path.join(details_dir, gridbox_files[0])
                print(f"\nFound gridbox file in ../details/ subdirectory: {gridbox_files[0]}")
            elif len(gridbox_files) > 1:
                print(f"\nFound multiple gridbox files in ../details/:")
                for i, f in enumerate(gridbox_files, 1):
                    print(f"  {i}. {f}")
                while True:
                    choice = input(f"Select file (1-{len(gridbox_files)}): ")
                    if choice.isdigit() and 1 <= int(choice) <= len(gridbox_files):
                        gridbox_file = os.path.join(details_dir, gridbox_files[int(choice) - 1])
                        break
                    print("Invalid selection.")
        
        # If not found automatically, prompt user
        while True:
            if gridbox_file is None:
                gridbox_file = input("Enter the path to the gridbox parameters file: ")
                if gridbox_file.startswith('"') and gridbox_file.endswith('"'):
                    gridbox_file = gridbox_file[1:-1]
            else:
                print(f"Reading gridbox parameters from: {gridbox_file}")
            
            if not os.path.exists(gridbox_file):
                print("That path does not appear to exist.")
                gridbox_file = None
                continue
                
            gridbox_params = read_gridbox_file(gridbox_file)
            if gridbox_params:
                break
            else:
                gridbox_file = None

    # At this point, gridbox_params should be populated
    if gridbox_params is None:
        print("\nError: Failed to obtain gridbox parameters. Exiting.")
        return

    # Display gridbox parameters
    print("\n" + "="*60)
    print("GRIDBOX PARAMETERS")
    print("="*60)
    print(f"Grid box size:   {gridbox_params['size_x']} × {gridbox_params['size_y']} × {gridbox_params['size_z']} Å")
    print(f"Grid box center: ({gridbox_params['center_x']}, {gridbox_params['center_y']}, {gridbox_params['center_z']})")

    # Get all PDBQT files in the directory
    pdbqt_files = [f for f in os.listdir(pwd) if f.endswith(".pdbqt") and not f.startswith(ligand_name)]

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
        
        # Check if this is a single structure or ensemble
        if 0 in ensemble_nums and len(ensemble_nums) == 1:
            single_count += 1
            print(f"\nSingle: {base_name}")
        else:
            # Count actual ensemble members (exclude 0)
            actual_ensemble_nums = [n for n in ensemble_nums if n != 0]
            if actual_ensemble_nums:
                ensemble_count += len(actual_ensemble_nums)
                print(f"\nEnsemble: {base_name}")
                print(f"  Members: {len(actual_ensemble_nums)} structures (ensemble_{min(actual_ensemble_nums)} to ensemble_{max(actual_ensemble_nums)})")
    
    print(f"\nTotal: {single_count} single structures, {ensemble_count} ensemble members")
    print("="*60)
    
    # Ask the user for the remaining configuration file information
    while True:        
        scoring = input("Enter the scoring function (ad4, vina [default] or vinardo): ")
        if not scoring:
            scoring = "vina"
        
        num_modes = input("Enter the number of modes (default: 5): ")
        if not num_modes:
            num_modes = "5"
        
        energy_range = input("Enter the energy range (default: 10): ")
        if not energy_range:
            energy_range = "10"
        
        exhaustiveness = input("Enter the exhaustiveness (default: 16; range = 8 - 32): ")
        if not exhaustiveness:
            exhaustiveness = "16"
        
        # Redisplay all of this information to the user and ask them to confirm that it is correct
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
            break
    
    # Ask about flexible residues
    while True:
        flex = input("\nAre there any flexible residues? (y/n): ").lower()
        while flex not in ["y", "n"]:
            print("Invalid input.")
            flex = input("Are there any flexible residues? (y/n): ").lower()
        
        flex_info = {}  # {(base_name, ensemble_num): has_flex_residues}
        
        if flex == "y":
            # Check for flex_residues.txt in ../details subdirectory first
            parent_dir = os.path.dirname(pwd)
            details_residues_file = os.path.join(parent_dir, "details", "flex_residues.txt")
            
            # Try automatic discovery
            if os.path.exists(details_residues_file):
                print(f"\nFound flex_residues.txt in ../details/ subdirectory")
                residues = details_residues_file
            elif os.path.exists(os.path.join(parent_dir, "flex_residues.txt")):
                print(f"\nFound flex_residues.txt in parent directory")
                residues = os.path.join(parent_dir, "flex_residues.txt")
            else:
                # Prompt user for the path
                residues = input("Enter the path to the flex_residues.txt file: ")
                if residues.startswith('"') and residues.endswith('"'):
                    residues = residues[1:-1]
                
                while not os.path.exists(residues):
                    residues = input(
                        "That path does not appear to exist.\nPlease enter the path to "
                        "the flex_residues.txt file: "
                    )
                    if residues.startswith('"') and residues.endswith('"'):
                        residues = residues[1:-1]
            
            print(f"Reading flexible residues from: {residues}")
            
            # Use the new parsing function
            try:
                flex_info = read_flexible_residues(residues)
                
                # Report what was loaded
                print(f"\nLoaded flexible residue data for {len(flex_info)} structure(s)")
                if flex_info:
                    structures_with_flex = sum(1 for has_flex in flex_info.values() if has_flex)
                    print(f"  Structures with flexible residues: {structures_with_flex}")
                    print(f"  Structures without flexible residues: {len(flex_info) - structures_with_flex}")
                    
            except Exception as e:
                print(f"\nError reading flexible residues file: {e}")
                print("Continuing without flexible residue information...")
                flex_info = {}
        
        # Categorize files
        flexible_configs = []
        rigid_configs = []
        
        for base_name in grouped_files:
            for ensemble_num in grouped_files[base_name]:
                file_info = grouped_files[base_name][ensemble_num]
                
                # Skip if this is the ligand file
                if base_name == ligand_name:
                    continue
                
                # Check if this structure has flexible residues
                # For ensemble structures, also check if base structure (ensemble_num=0) has flex info
                has_flex = flex_info.get((base_name, ensemble_num), False)
                if not has_flex and ensemble_num != 0:
                    # Check if base structure has flex residues (applies to all ensemble members)
                    has_flex = flex_info.get((base_name, 0), False)
                
                if has_flex:
                    # Must have both rigid and flex files
                    if 'rigid' in file_info and 'flex' in file_info:
                        flexible_configs.append((base_name, ensemble_num, file_info))
                    else:
                        if ensemble_num == 0:
                            print(f"Warning: {base_name} marked as flexible but missing rigid/flex files")
                        else:
                            print(f"Warning: {base_name}_ensemble_{ensemble_num} marked as flexible but missing rigid/flex files")
                else:
                    # Use full file if available, otherwise rigid
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
                if ens_num == 0:
                    print(f"  {base_name}")
                else:
                    print(f"  {base_name}_ensemble_{ens_num}")
        else:
            print("  (none)")
        
        print("\nRIGID RECEPTORS:")
        if rigid_configs:
            for base_name, ens_num, _ in rigid_configs:
                if ens_num == 0:
                    print(f"  {base_name}")
                else:
                    print(f"  {base_name}_ensemble_{ens_num}")
        else:
            print("  (none)")
        print("="*60)
        
        confirm = input("\nIs this categorization correct? (y/n): ").lower()
        while confirm not in ["y", "n"]:
            print("Invalid input.")
            confirm = input("Is this categorization correct? (y/n): ").lower()
        
        if confirm == "y":
            break
    
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
    print(f"Files saved to: {os.path.abspath(pwd)}")


if __name__ == "__main__":
    get_config_files()