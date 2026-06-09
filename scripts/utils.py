"""
Shared utility functions for the cross-species docking pipeline.
"""
import os
import json
import shutil
import numpy as np
from scipy.optimize import linear_sum_assignment

# ---------------------------------------------------------------------------
# Project config helpers
# ---------------------------------------------------------------------------

CONFIG_FILENAME = "project_config.json"


def find_project_dir(start_path=None):
    """
    Walk up from start_path looking for project_config.json.
    Once found, reads the 'project_dir' key stored inside it so the answer
    is always authoritative regardless of which subdirectory a script is run
    from.  Falls back to the directory that contains the file if the key is
    missing or the stored path no longer exists.
    Returns the absolute project root path, or None if no config file is found.
    """
    search = os.path.abspath(start_path if start_path else os.getcwd())
    for _ in range(10):
        cfg_path = os.path.join(search, CONFIG_FILENAME)
        if os.path.isfile(cfg_path):
            try:
                with open(cfg_path) as f:
                    cfg = json.load(f)
                stored = cfg.get("project_dir")
                if stored and os.path.isdir(stored):
                    return stored
            except (json.JSONDecodeError, OSError):
                pass
            return search  # config found but unreadable or key missing
        parent = os.path.dirname(search)
        if parent == search:
            break
        search = parent
    return None


def load_config(project_dir):
    """Load project_config.json from project_dir. Returns empty dict if absent."""
    path = os.path.join(project_dir, CONFIG_FILENAME)
    if os.path.isfile(path):
        with open(path) as f:
            return json.load(f)
    return {}


def save_config(project_dir, config):
    """Write config dict to project_config.json in project_dir."""
    with open(os.path.join(project_dir, CONFIG_FILENAME), "w") as f:
        json.dump(config, f, indent=2)


def get_project_paths(project_dir):
    """Return a dict of all standard subdirectory paths for a project."""
    return {
        "project_dir":         project_dir,
        "prepared_structures": os.path.join(project_dir, "prepared_structures"),
        "pdbqt_files":         os.path.join(project_dir, "pdbqt_files"),
        "docking_results":     os.path.join(project_dir, "docking_results"),
        "models":              os.path.join(project_dir, "docking_results", "models"),
        "results":             os.path.join(project_dir, "results"),
    }


def resolve_project_dir(hint=None):
    """
    Locate the project root directory and return its absolute path.

    Walks up from hint (or CWD) looking for project_config.json, then reads
    the 'project_dir' value stored inside it.  Confirms the detected path with
    the user before using it.  Falls back to a manual prompt if no config file
    is found or the user rejects the detected path.
    """
    found = find_project_dir(hint)
    if found:
        print(f"\nFound project directory: {found}")
        ans = input("Use this? (y/n): ").lower().strip()
        if ans == "y":
            return found

    while True:
        path = input("Enter the project directory path: ").strip().strip('"')
        if os.path.isdir(path):
            return os.path.abspath(path)
        print("Directory not found. Please try again.")


def resolve_reference_pdb(config, project_dir):
    """
    Locate the reference PDB file from the project config.

    Searches first in prepared_structures/ (post-prep version with ligand
    preserved), then in project_dir/ (original input file).  Returns the full
    absolute path, or None if the config has no reference_pdb entry or the
    file cannot be found.
    """
    basename = config.get("reference_pdb")
    if not basename:
        return None
    paths = get_project_paths(project_dir)
    for search_dir in [paths["prepared_structures"], project_dir]:
        candidate = os.path.join(search_dir, basename)
        if os.path.exists(candidate):
            return candidate
    return None


def euclidean3d(v1, v2):
    """Calculate Euclidean distance between two 3D points."""
    if not (len(v1) == 3 and len(v2) == 3):
        return None
    return np.sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def centroid(coords):
    """Calculate the centroid from a 3D point cloud."""
    return list(
        map(
            np.mean,
            ([c[0] for c in coords], [c[1] for c in coords], [c[2] for c in coords]),
        )
    )


def parse_hetatm_coords(lines, lig_name):
    """
    Extract heavy-atom coordinates and element symbols from HETATM records.

    Works identically on both PDB and PDBQT files because both formats share
    the same fixed-column layout for the coordinate block (columns 31–54).
    Hydrogen atoms (element symbol 'H') are skipped.

    Parameters
    ----------
    lines : list[str]
        Raw text lines from a PDB or PDBQT file (or a MODEL block subset).
    lig_name : str
        Three-letter residue name of the ligand to extract (e.g. 'DHT').

    Returns
    -------
    coords : np.ndarray, shape (N, 3)
        Heavy-atom x/y/z coordinates.
    elements : list[str]
        One-letter element symbol per heavy atom.
    """
    coords, elements = [], []
    for line in lines:
        if line.startswith("HETATM") and line[17:20].strip() == lig_name:
            atom_name = line[12:16].strip()
            element = atom_name[:1]
            if element == "H":
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            elements.append(element)
            coords.append([x, y, z])
    return np.array(coords) if coords else np.empty((0, 3)), elements


def rmsd_hungarian(ref_coords, query_coords):
    """
    Compute RMSD between two coordinate arrays using the Hungarian algorithm
    for optimal atom assignment.

    Using sorted-axis matching (np.sort) destroys x/y/z triplet relationships
    and produces geometrically meaningless RMSDs for symmetric or reordered
    ligands.  The Hungarian algorithm finds the minimum-cost bipartite
    matching between atoms, guaranteeing the correct geometric RMSD regardless
    of atom ordering in the input files.

    Parameters
    ----------
    ref_coords : np.ndarray, shape (N, 3)
        Reference heavy-atom coordinates.
    query_coords : np.ndarray, shape (N, 3)
        Query heavy-atom coordinates (same number of atoms as ref).

    Returns
    -------
    float
        RMSD in Ångströms after optimal atom matching.
    """
    cost = np.linalg.norm(ref_coords[:, None, :] - query_coords[None, :, :], axis=2)
    row_ind, col_ind = linear_sum_assignment(cost)
    matched = query_coords[col_ind]
    diff = ref_coords[row_ind] - matched
    return float(np.sqrt(np.sum(diff ** 2) / len(row_ind)))


def check_tools(tools):
    """
    Verify that external command-line tools are available on PATH.

    Parameters
    ----------
    tools : list[str]
        Tool names to check (e.g. ["muscle", "vina", "plip"]).

    Raises
    ------
    SystemExit
        If any required tool is missing, prints the missing names and exits.
    """
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        print("ERROR: The following required tools were not found on PATH:")
        for t in missing:
            print(f"  - {t}")
        print("\nPlease install the missing tools and ensure they are accessible from PATH.")
        raise SystemExit(1)
