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
    Read flexible residues from the summary file created by mds_structure_prep.py
    
    The expected format is:
    FLEXIBLE RESIDUES SUMMARY
    ==========================================================
    
    structure_name: [resid1,resid2,resid3]
    another_structure: [resid4,resid5]
    
    ==========================================================
    Total structures analyzed: N
    Structures with flexible residues: M
    
    Parameters:
    -----------
    residues_file : str
        Path to the flex_residues.txt file
    
    Returns:
    --------
    dict: Dictionary with structure base names as keys and flexible residues lists as values
    """
    flex_residues = {}
    
    with open(residues_file, "r") as f:
        lines = f.readlines()
    
    # Skip header lines and statistics footer
    for line in lines:
        line = line.strip()
        
        # Skip empty lines, header, separator lines, and statistics
        if not line or line.startswith('=') or line.startswith('FLEXIBLE') or \
           line.startswith('Total') or line.startswith('Structures with'):
            continue
        
        # Parse structure_name: [resid1,resid2,resid3] format
        if ':' in line and '[' in line and ']' in line:
            try:
                structure_name, residues_str = line.split(':', 1)
                structure_name = structure_name.strip()
                
                # Extract residues from brackets
                residues_str = residues_str.strip()
                residues_str = residues_str.split('[')[1].split(']')[0].strip()
                
                # Store as list of residue IDs (empty list if no residues)
                if residues_str:
                    residues = [r.strip() for r in residues_str.split(',')]
                else:
                    residues = []
                
                flex_residues[structure_name] = residues
                
            except (ValueError, IndexError) as e:
                print(f"Warning: Could not parse line: {line}")
                print(f"  Error: {e}")
                continue

    # List all species that had flexible residue data loaded
    print(f"\nFlexible residue data loaded for {len(flex_residues)} structure(s):")
    for structure in flex_residues:
        print(f"  {structure}: {len(flex_residues[structure])} flexible residue(s)")
    
    return flex_residues


def get_base_structure_name(pdb_filename):
    """
    Extract the base structure name from a PDB filename.
    Handles ensemble structures by removing the _ensemble_N suffix.
    
    Examples:
    ---------
    'structure.pdb' -> 'structure'
    'structure_ensemble_1.pdb' -> 'structure'
    'structure_ensemble_10.pdb' -> 'structure'
    'my_structure_fixed.pdb' -> 'my_structure_fixed'
    'my_structure_fixed_ensemble_3.pdb' -> 'my_structure_fixed'
    
    Parameters:
    -----------
    pdb_filename : str
        PDB filename (with or without path)
    
    Returns:
    --------
    str: Base structure name without .pdb extension and without _ensemble_N suffix
    """
    # Remove .pdb extension
    base_name = os.path.basename(pdb_filename).replace('.pdb', '')
    
    # Check if this is an ensemble structure
    if '_ensemble_' in base_name:
        # Split on _ensemble_ and take everything before it
        parts = base_name.split('_ensemble_')
        base_name = parts[0]
    
    return base_name


def format_flexible_residues(pdb_file, residue_ids):
    """
    Format flexible residues with chain identifiers from the PDB file.
    
    This reads the PDB file to extract chain information for each residue ID,
    then formats them as required by mk_prepare_receptor (chain:resid format).
    
    Parameters:
    -----------
    pdb_file : str
        Path to the PDB file
    residue_ids : list
        List of residue IDs as strings (e.g., ['123', '456', '789'])
    
    Returns:
    --------
    str: Formatted flexible residues string (e.g., 'A:123,A:456,B:789')
         Returns empty string if no residues or if none could be formatted
    """
    if not residue_ids:
        return ""
    
    # Read PDB file to get chain information
    formatted_residues = []
    residues_to_find = set(residue_ids.copy())
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    resid = line[22:26].strip()
                    if resid in residues_to_find:
                        chain = line[21]
                        formatted_res = f"{chain}:{resid}"
                        if formatted_res not in formatted_residues:
                            formatted_residues.append(formatted_res)
                        residues_to_find.discard(resid)
                        
                        # Stop if we've found all residues
                        if not residues_to_find:
                            break
    
    except Exception as e:
        print(f"Warning: Error reading PDB file {pdb_file}: {e}")
        return ""
    
    if residues_to_find:
        print(f"Warning: Could not find chain info for residues: {sorted(residues_to_find)}")
    
    return ','.join(formatted_residues)


def main():
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from utils import resolve_project_dir, get_project_paths

    project_dir = resolve_project_dir()
    paths = get_project_paths(project_dir)
    pdb_dir = paths["prepared_structures"]

    if not os.path.isdir(pdb_dir):
        print(f"Error: prepared_structures/ not found at {pdb_dir}")
        print("Please run mds_structure_prep.py first.")
        return

    os.chdir(pdb_dir)
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith(".pdb")]
    print(f"\nFound {len(pdb_files)} PDB file(s)")

    # ── Flexible residues (auto-detected from results/flex_residues.txt) ──────
    flex_residues_dict = {}
    details_residues_file = os.path.join(paths["results"], "flex_residues.txt")

    if os.path.exists(details_residues_file):
        print(f"Found flex_residues.txt in results/ — loading flexible residue data...")
        flex_residues_dict = get_flexible_residues(details_residues_file)

        structure_names = set(get_base_structure_name(pdb) for pdb in pdb_files)
        missing_structures = structure_names - set(flex_residues_dict.keys())
        if missing_structures:
            print(f"\nWarning: Flexible residue data is missing for:")
            for struct in missing_structures:
                print(f"  {struct}")
            print("These structures will be prepared as rigid receptors.")
    else:
        print("No flex_residues.txt found in results/ — all receptors will be prepared as rigid.")

    # Process all PDB files
    processed_count = 0
    skipped_count = 0
    
    for pdb in pdb_files:
        if pdb.startswith("ref_"):
            print(f"\nSkipping reference structure {pdb}...")
            skipped_count += 1
            continue
        
        # Get base structure name (handles ensemble structures)
        base_name = get_base_structure_name(pdb)
        pdb_name = pdb.replace(".pdb", "")
        
        print(f"\nPreparing {pdb_name}...")
        if '_ensemble_' in pdb:
            print(f"  (Ensemble structure - using flexible residues from base structure: {base_name})")

        # Check if we have flexible residues for this structure's base name
        if base_name in flex_residues_dict:
            residue_ids = flex_residues_dict[base_name]
            
            if residue_ids:
                # Format residues with chain information
                formatted_residues = format_flexible_residues(pdb, residue_ids)
                
                if formatted_residues:
                    print(f"  Applying {len(residue_ids)} flexible residue(s): {formatted_residues}")
                    prepare_receptor(pdb, pdb_name, formatted_residues)
                else:
                    print(f"  Warning: Could not format flexible residues. Preparing rigid receptor...")
                    prepare_receptor(pdb, pdb_name)
            else:
                print(f"  No flexible residues specified for {base_name}. Preparing rigid receptor...")
                prepare_receptor(pdb, pdb_name)
        else:
            # No flexible residue data for this structure
            if flex_residues_dict:  # Only warn if we loaded flex data but this structure isn't in it
                print(f"  No flexible residue data found for {base_name}. Preparing rigid receptor...")
            prepare_receptor(pdb, pdb_name)
        
        processed_count += 1

    # Move PDBQT files to output directory (sibling of prepared_structures/)
    pdbqt_out = os.path.normpath(os.path.join(pdb_dir, "..", "pdbqt_files"))
    os.makedirs(pdbqt_out, exist_ok=True)

    pdbqt_count = 0
    for f in os.listdir(pdb_dir):
        if f.endswith(".pdbqt"):
            os.replace(f, os.path.join(pdbqt_out, f))
            pdbqt_count += 1

    # Print summary
    print("\n" + "="*60)
    print("Processing Complete")
    print("="*60)
    print(f"\nStructures processed: {processed_count}")
    print(f"Structures skipped (references): {skipped_count}")
    print(f"PDBQT files created: {pdbqt_count}")
    print(f"\nSuccessfully prepared structures have been moved to: {pdbqt_out}")

if __name__ == "__main__":
    main()