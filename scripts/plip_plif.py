"""
Integrated PLIP Analysis and PLIF Generation Pipeline

Workflow:
1. Runs PLIP on all PDB files in specified directory
2. Processes protonated structures with OpenBabel
3. Organizes PLIP outputs (XML and protonated PDB files)
4. Automatically generates PLIFs for all test structures vs reference
5. Calculates Tanimoto similarity coefficients
6. Creates organized output with individual comparisons and summary table
"""

import os
import sys
import shutil
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import euclidean3d, check_tools


def run_plip_analysis(dir_path, results_dir=None):
    """
    Run PLIP on all PDB files in the specified directory and organize results.
    
    Parameters:
    -----------
    dir_path : str
        Path to directory containing PDB files to analyze
    
    Returns:
    --------
    str
        Path to the PLIP_results directory containing XML and protonated PDB files
    
    Methodology:
    -----------
    PLIP (Protein-Ligand Interaction Profiler) detects non-covalent interactions
    including hydrogen bonds, hydrophobic contacts, pi-stacking, salt bridges, etc.
    The -xv flag generates XML output with detailed interaction geometry, and --name
    ensures consistent file naming for downstream processing.
    
    OpenBabel is used to add hydrogens to the protonated structures, which is
    necessary for accurate representation of hydrogen bonding and protonation states
    that may vary between species or pH conditions.
    """
    print("\n" + "=" * 70)
    print("Running PLIP Analysis")
    print("=" * 70)
    
    # Change to specified directory
    os.chdir(dir_path)
    
    # Run PLIP on all PDB files
    pdb_files = [f for f in os.listdir(dir_path) if f.endswith(".pdb")]
    print(f"\nFound {len(pdb_files)} PDB files to analyze")
    
    for file in pdb_files:
        base_name = os.path.splitext(file)[0]
        try:
            print(f"  Running PLIP on {file}...")
            subprocess.run(
                ["plip", "-f", file, "-xv", "--name", base_name],
                check=True,
                capture_output=True
            )
        except subprocess.CalledProcessError as e:
            print(f"  Error processing {file}: {e}")
        except Exception as e:
            print(f"  Unexpected error processing {file}: {e}")
    
    # Process protonated files with OpenBabel
    print("\nProcessing protonated structures with OpenBabel...")
    protonated_files = [f for f in os.listdir(dir_path) if f.endswith('_protonated.pdb')]
    for file in protonated_files:
        print(f"  Adding hydrogens to {file}...")
        subprocess.run(
            ["obabel", file, "-o", "pdb", "-O", file, "-h"],
            capture_output=True
        )
    
    # Create and organize results directory (caller can override via results_dir)
    if results_dir is None:
        results_dir = os.path.join(dir_path, "PLIP_results")
    os.makedirs(results_dir, exist_ok=True)
    
    print("\nOrganizing files into PLIP_results directory...")
    
    # Move XML and protonated PDB files
    files_moved = 0
    for file in os.listdir(dir_path):
        if file.endswith(".xml") or "_protonated.pdb" in file:
            src = os.path.join(dir_path, file)
            dst = os.path.join(results_dir, file)
            try:
                shutil.move(src, dst)
                files_moved += 1
            except Exception as e:
                print(f"  Warning: Could not move {file}: {e}")
    
    print(f"Moved {files_moved} files to PLIP_results folder")
    
    return results_dir


def parse_plip_xml(xml_file):
    """
    Parse PLIP XML output and return DataFrame with residue number, 
    residue type, and interaction type for each interaction.
    
    Methodology:
    -----------
    PLIP organizes interactions by type in the XML structure. We extract all
    interaction types and represent them at the residue level (not atom level)
    because residue-level patterns are more stable across species than specific
    atom positions, which can vary due to side-chain conformations.
    
    Duplicates are removed because we want to know IF a residue participates
    in each interaction type, not HOW MANY times. This creates a binary
    fingerprint that's more robust to minor structural differences that can
    exist between species.
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    data = []
    
    # Try both capitalisation variants emitted by different PLIP versions,
    # then fall back to scanning all binding sites.
    # Explicit `is not None` checks are required — Element truth-testing is deprecated.
    binding_site = root.find('bindingsite[@id="1"][@has_interactions="True"]')
    if binding_site is None:
        binding_site = root.find('bindingsite[@id="1"][@has_interactions="true"]')
    if binding_site is None:
        binding_site = next(
            (bs for bs in root.findall("bindingsite")
             if bs.get("has_interactions", "").lower() == "true"),
            None,
        )
    if binding_site is not None:
        interactions = binding_site.find("interactions")
        
        for interaction_type in [
            "hydrophobic_interactions",
            "hydrogen_bonds",
            "water_bridges",
            "salt_bridges",
            "pi-stacks",
            "pi_cation_interactions",
            "halogen_bonds",
            "metal_complexes",
        ]:
            interaction_elem = interactions.find(interaction_type)
            if interaction_elem is not None:
                for interaction in interaction_elem:
                    resnr = interaction.find("resnr")
                    if resnr is not None:
                        resnr = resnr.text
                        restype = interaction.find("restype").text
                        data.append({
                            "resnr": resnr,
                            "restype": restype,
                            "interaction_type": interaction_type,
                        })
    
    df = pd.DataFrame(data)
    # Remove duplicates to count interactions on a per residue basis
    df = df.drop_duplicates(keep="first")
    return df


def get_vdw_contacts(pdb_file, ligand_name):
    """
    Calculate van der Waals contacts between protein and ligand.
    
    Parameters:
    -----------
    pdb_file : str
        Path to protonated PDB file
    ligand_name : str
        Ligand identifier as found in HETATM records
    
    Returns:
    --------
    DataFrame
        van der Waals contacts with residue number, type, and interaction type
    
    Methodology:
    -----------
    VDW contacts are detected by comparing inter-atomic distances to the sum of
    VDW radii plus a 0.6 Å tolerance. This captures weak contacts that may not
    form classical interaction types but still contribute to binding affinity.
    
    This is particularly important for cross-species comparisons because even if
    specific H-bonds or salt bridges differ due to sequence variation, conserved
    VDW contacts can indicate similar binding modes. The tolerance of 0.6 Å
    accounts for thermal fluctuations and coordinate uncertainty.
    """
    vdw_radii = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}
    protein_coordinates = []
    ligand_coordinates = []
    
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                atom = line[12:16].strip()[0]
                resnr = line[22:26].strip()
                restype = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                protein_coordinates.append([atom, resnr, restype, x, y, z])
            elif line.startswith("HETATM") and line[17:20].strip() == ligand_name:
                atom = line[12:16].strip()[0]
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ligand_coordinates.append([atom, x, y, z])
    
    if protein_coordinates and ligand_coordinates:
        all_rows = []
        for lig_coord in ligand_coordinates:
            for prot_coord in protein_coordinates:
                dist = euclidean3d(
                    [lig_coord[1], lig_coord[2], lig_coord[3]],
                    [prot_coord[3], prot_coord[4], prot_coord[5]]
                )
                all_rows.append({
                    "HETATM": lig_coord[0],
                    "LIG x": lig_coord[1],
                    "LIG y": lig_coord[2],
                    "LIG z": lig_coord[3],
                    "ATOM": prot_coord[0],
                    "resnr": prot_coord[1],
                    "restype": prot_coord[2],
                    "PROT x": prot_coord[3],
                    "PROT y": prot_coord[4],
                    "PROT z": prot_coord[5],
                    "DIST": dist
                })
        
        df = pd.DataFrame(all_rows)
        
        # Calculate vdw_radii and interactions
        df["vdw_radii"] = df.apply(
            lambda row: vdw_radii.get(row["HETATM"][0], 0) + 
                       vdw_radii.get(row["ATOM"][0], 0) + 0.6,
            axis=1
        )
        df["vdw_interaction"] = df["DIST"] < df["vdw_radii"]
        
        vdw_df = df[df["vdw_interaction"]].copy()
        vdw_df = vdw_df.groupby(["resnr", "restype"]).first().reset_index()
        vdw_df = vdw_df[["resnr", "restype"]].assign(interaction_type="vdw_contact")
    else:
        vdw_df = pd.DataFrame(columns=["resnr", "restype", "interaction_type"])
    
    return vdw_df


def merge_plifs(ref_df, test_df):
    """
    Merge reference and test PLIFs to create comparison fingerprint.
    
    Parameters:
    -----------
    ref_df : DataFrame
        Reference PLIF with resnr, restype, interaction_type
    test_df : DataFrame
        Test PLIF with resnr, restype, interaction_type
    
    Returns:
    --------
    DataFrame
        Merged fingerprint with binary columns indicating presence in ref/test
    
    Methodology:
    -----------
    Matching key: restype + interaction_type (NOT resnr).

    In a cross-species pipeline each receptor has its own native residue
    numbering, so a conserved binding residue (e.g. LEU in the androgen
    receptor binding pocket) will have a different sequence number in every
    species.  Merging on residue number would therefore never find a match and
    would always produce Tanimoto = 0.

    By collapsing to unique (restype, interaction_type) pairs before merging,
    we ask the biologically meaningful question: "does this species form a
    hydrogen bond with a THR residue?" rather than "does residue 710 form a
    hydrogen bond?".  This is the standard approach for cross-species PLIF
    comparison.

    Rows in the output retain the reference residue number (resnr) for display
    purposes so the output file can still be used to trace interactions back to
    specific positions in the reference structure.
    """
    # Collapse to one row per (restype, interaction_type) per structure.
    # If multiple same-type residues make the same interaction class, only the
    # first occurrence is kept — each feature is counted once (binary).
    ref_uniq = (
        ref_df
        .drop_duplicates(subset=["restype", "interaction_type"], keep="first")
        .rename(columns={
            "resnr": "resnr_ref",
            "restype": "restype_ref",
            "interaction_type": "interaction_type_ref",
        })
    )

    test_uniq = (
        test_df
        .drop_duplicates(subset=["restype", "interaction_type"], keep="first")
        .rename(columns={
            "resnr": "resnr_test",
            "restype": "restype_test",
            "interaction_type": "interaction_type_test",
        })
    )

    merged_df = pd.merge(
        ref_uniq,
        test_uniq,
        left_on=["restype_ref", "interaction_type_ref"],
        right_on=["restype_test", "interaction_type_test"],
        how="outer",
    )

    # Binary presence flags — set before any fillna so NaN correctly signals absence
    merged_df = merged_df.assign(
        ref=merged_df["resnr_ref"].notnull().astype(int),
        test=merged_df["resnr_test"].notnull().astype(int),
    )

    # Fill display columns from whichever side has data
    merged_df["resnr_ref"] = merged_df["resnr_ref"].fillna(merged_df["resnr_test"])
    merged_df["restype_ref"] = merged_df["restype_ref"].fillna(merged_df["restype_test"])
    merged_df["interaction_type_ref"] = merged_df["interaction_type_ref"].fillna(
        merged_df["interaction_type_test"]
    )

    result_df = merged_df[
        ["resnr_ref", "restype_ref", "interaction_type_ref", "ref", "test"]
    ].rename(columns={
        "resnr_ref": "resnr",
        "restype_ref": "restype",
        "interaction_type_ref": "interaction",
    })

    return result_df


def calculate_tanimoto(plif_df):
    """
    Calculate Tanimoto similarity coefficient from merged PLIF DataFrame.
    
    Methodology:
    -----------
    Tanimoto coefficient (also called Jaccard index for binary data) measures
    similarity as: (shared features) / (total unique features)

        Tanimoto = |A ∩ B| / |A ∪ B|

    The coefficient ranges from 0 (no overlap) to 1 (perfect match), providing
    an intuitive measure of binding site similarity across species.
    
    """
    both  = int(((plif_df["ref"] == 1) & (plif_df["test"] == 1)).sum())
    either = int(((plif_df["ref"] == 1) | (plif_df["test"] == 1)).sum())
    tanimoto = both / either if either > 0 else 0.0
    return round(tanimoto, 3)


def generate_heatmap(all_plifs, ref_df, ref_name, output_dir, tanimoto_scores=None):
    """
    Generate a binary PLIF heatmap: species on the y-axis, reference binding
    pocket interactions on the x-axis (top).

    Parameters:
    -----------
    all_plifs : dict
        {pdb_name: merged_plif_df} — one entry per species (best model already
        selected by the caller)
    ref_df : DataFrame
        Reference PLIF (columns: resnr, restype, interaction_type)
    ref_name : str
        Display name for the reference structure
    output_dir : str
        Directory in which to save plif_heatmap.png
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.colors import ListedColormap
    except ImportError:
        print("  Warning: matplotlib not available — skipping heatmap.")
        return

    # Keys match the raw PLIP/VDW interaction_type strings in the data
    itype_colors = {
        "hydrophobic_interactions": "#df7111",
        "hydrogen_bonds":           "#2980b9",
        "water_bridges":            "#00bcd4",
        "salt_bridges":             "#e74c3c",
        "pi-stacks":                "#9b59b6",
        "pi_cation_interactions":   "#f1c40f",
        "halogen_bonds":            "#27ae60",
        "metal_complexes":          "#789da0",
        "vdw_contact":              "#455152",
    }

    # Human-readable labels for the legend
    itype_labels = {
        "hydrophobic_interactions": "Hydrophobic",
        "hydrogen_bonds":           "Hydrogen bond",
        "water_bridges":            "Water bridge",
        "salt_bridges":             "Salt bridge",
        "pi-stacks":                "Pi-stacking",
        "pi_cation_interactions":   "Pi-cation interaction",
        "halogen_bonds":            "Halogen bond",
        "metal_complexes":          "Metal complex",
        "vdw_contact":              "van der Waals",
    }

    # ── Columns: grouped by interaction type, then sorted by residue number ───
    itype_order = {k: i for i, k in enumerate(itype_colors)}
    ref_unique = ref_df.drop_duplicates(subset=["resnr", "restype", "interaction_type"])
    try:
        ref_sorted = ref_unique.copy()
        ref_sorted["_itype_order"] = ref_sorted["interaction_type"].map(itype_order).fillna(99)
        ref_sorted["_resnr_int"] = ref_sorted["resnr"].astype(int)
        ref_sorted = ref_sorted.sort_values(
            ["_itype_order", "_resnr_int"]
        ).drop(columns=["_itype_order", "_resnr_int"]).reset_index(drop=True)
    except (ValueError, TypeError):
        ref_sorted = ref_unique.copy()
        ref_sorted["_itype_order"] = ref_sorted["interaction_type"].map(itype_order).fillna(99)
        ref_sorted = ref_sorted.sort_values(
            ["_itype_order", "resnr"]
        ).drop(columns="_itype_order").reset_index(drop=True)

    col_keys = list(
        zip(ref_sorted["resnr"], ref_sorted["restype"], ref_sorted["interaction_type"])
    )
    col_labels = [
        f"{resnr}{restype}"
        for resnr, restype, _ in col_keys
    ]

    # ── Rows: reference first, then species sorted by Tanimoto (high → low) ───
    sorted_pdbs = sorted(
        all_plifs.keys(),
        key=lambda p: (tanimoto_scores or {}).get(p, 0),
        reverse=True,
    )
    species_labels = [pdb.split("_")[0] for pdb in sorted_pdbs]
    row_labels = ["Reference"] + species_labels

    matrix = [np.ones(len(col_keys), dtype=int)]  # reference row — always all present
    for pdb in sorted_pdbs:
        plif_df = all_plifs[pdb]
        row = []
        for _, restype, itype in col_keys:
            mask = (plif_df["restype"] == restype) & (plif_df["interaction"] == itype)
            row.append(int(plif_df.loc[mask, "test"].iloc[0]) if mask.any() else 0)
        matrix.append(row)

    matrix = np.array(matrix)
    n_rows, n_cols = matrix.shape

    # ── Figure sizing ─────────────────────────────────────────────────────────
    fig_w = max(14, n_cols * 0.7 + 5)
    fig_h = max(6,  n_rows * 0.7 + 3)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = ListedColormap(["#D81B60", "#73B2E9"])  # pink = absent, blue = present
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1)

    # X-axis on bottom
    ax.xaxis.set_label_position("bottom")
    ax.xaxis.tick_bottom()
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, rotation=90, ha="center", va="top", fontsize=16)
    ax.set_xlabel("Binding pocket residue", fontsize=16, labelpad=10)

    # Color each x-axis tick label by its interaction type
    fig.canvas.draw()
    for tick, (_, _, itype) in zip(ax.get_xticklabels(), col_keys):
        tick.set_color(itype_colors.get(itype, "black"))

    # Y-axis (left)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=16)
    ax.set_ylabel("Species", fontsize=16, labelpad=10)

    # Grid
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", size=0)

    # ── Colored column-group borders by interaction type ──────────────────────
    # Identify contiguous runs of the same interaction type
    col_groups = []
    if col_keys:
        cur_itype = col_keys[0][2]
        cur_start = 0
        for i, (_, _, itype) in enumerate(col_keys[1:], 1):
            if itype != cur_itype:
                col_groups.append((cur_start, i - 1, cur_itype))
                cur_itype = itype
                cur_start = i
        col_groups.append((cur_start, len(col_keys) - 1, cur_itype))

    for grp_start, grp_end, grp_itype in col_groups:
        ax.add_patch(mpatches.Rectangle(
            (grp_start - 0.5, -0.5),
            grp_end - grp_start + 1,
            n_rows,
            linewidth=5,
            edgecolor=itype_colors.get(grp_itype, "black"),
            facecolor="none",
            zorder=4,
            clip_on=False,
        ))

    # Separator after reference row
    ax.axhline(y=0.5, color="black", linewidth=2)

    # ── Right y-axis: Tanimoto scores ─────────────────────────────────────────
    if tanimoto_scores:
        tanimoto_labels = ["—"] + [
            f"{tanimoto_scores.get(pdb, ''):.3f}" for pdb in sorted_pdbs
        ]
        ax2 = ax.twinx()
        ax2.set_ylim(ax.get_ylim())
        ax2.set_yticks(range(n_rows))
        ax2.set_yticklabels(tanimoto_labels, fontsize=16)
        ax2.set_ylabel("Tanimoto similarity", fontsize=16, labelpad=10)
        ax2.tick_params(axis="y", length=0)

    # ── Legends (stacked, shifted right to clear the Tanimoto axis) ──────────
    # Legend 1: presence / absence
    seen_itypes = dict.fromkeys(itype for _, _, itype in col_keys)
    itype_handles = [
        mpatches.Patch(
            color=itype_colors.get(itype, "black"),
            label=itype_labels.get(itype, itype.replace("_", " ")),
        )
        for itype in seen_itypes
    ]
    presence_legend = ax.legend(
        handles=[
            mpatches.Patch(color="#73B2E9", label="Present"),
            mpatches.Patch(color="#D81B60", label="Absent"),
        ],
        loc="upper left",
        bbox_to_anchor=(1.06, 1.0),
        frameon=False,
        fontsize=14,
        title="Interaction",
        title_fontsize=16,
    )
    ax.add_artist(presence_legend)

    # Legend 2: interaction type colors
    ax.legend(
        handles=itype_handles,
        loc="upper left",
        bbox_to_anchor=(1.06, 0.78),
        frameon=False,
        fontsize=14,
        title="Interaction type",
        title_fontsize=16,
    )

    plt.tight_layout()
    out_path = os.path.join(output_dir, "plif_heatmap.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nHeatmap saved to: {out_path}")


def generate_plifs(plip_results_dir, ligand_name, ref_pdb=None, summary_dir=None):
    """
    Generate PLIFs for all test structures compared to reference.
    
    Parameters:
    -----------
    plip_results_dir : str
        Path to directory containing PLIP XML and protonated PDB files
    ligand_name : str
        Ligand identifier as found in HETATM records
    ref_pdb : str, optional
        Name of reference PDB file (if None, will prompt user)
    
    Returns:
    --------
    DataFrame
        Summary table with Tanimoto similarities for all comparisons
    
    Methodology:
    -----------
    This function implements the core PLIF comparison workflow.
    
    Individual PLIF comparison files are saved for further inspection.
    """
    print("\n" + "=" * 70)
    print("Generating Protein-Ligand Interaction Fingerprints")
    print("=" * 70)
    
    os.chdir(plip_results_dir)
    
    # Find reference PDB file if not specified
    if ref_pdb is None:
        pdb_files = [i for i in os.listdir(plip_results_dir) 
                     if i.endswith("_protonated.pdb")]
        
        # Auto-detect reference file by "ref_" prefix
        ref_candidates = [f for f in pdb_files if f.startswith("ref_")]

        if len(ref_candidates) == 1:
            ref_pdb = ref_candidates[0]
            print(f"\nAuto-detected reference PDB: {ref_pdb}")
        elif len(ref_candidates) > 1:
            print(f"\nMultiple ref_ files found: {ref_candidates}")

        if ref_pdb is None:
            print("\nAvailable protonated PDB files:")
            for i, file in enumerate(pdb_files, 1):
                print(f"  {i}. {file}")
            ref_pdb = input("\nEnter the name of the protonated reference PDB file: ")
            if not ref_pdb.endswith(".pdb"):
                ref_pdb += ".pdb"
            while not os.path.exists(ref_pdb):
                print('Error: This file does not appear to exist.')
                ref_pdb = input("Enter the file name of the protonated reference PDB file: ")
                if not ref_pdb.endswith(".pdb"):
                    ref_pdb += ".pdb"
    
    # Get reference base name
    ref_pdb_name = ref_pdb.replace("_protonated.pdb", "")
    ref_xml = ref_pdb_name + ".xml"
    
    print(f"\nProcessing reference structure: {ref_pdb_name}")
    print("=" * 70)
    
    # Parse reference PLIP data
    print(f"Parsing PLIP XML report for {ref_pdb_name}...")
    ref_plip_df = parse_plip_xml(ref_xml)
    
    # Get reference van der Waals contacts
    print(f"Calculating van der Waals contacts for {ref_pdb_name}...")
    ref_vdw_df = get_vdw_contacts(ref_pdb, ligand_name)
    
    # Combine reference PLIP and VDW data
    ref_df = pd.concat([ref_plip_df, ref_vdw_df]).reset_index(drop=True)
    
    print(f"\nReference PLIF contains {len(ref_df)} interactions")
    print(f"  - PLIP-detected: {len(ref_plip_df)}")
    print(f"  - VDW contacts: {len(ref_vdw_df)}")
    
    print("\nProcessing test structures...")
    print("=" * 70)

    # Process all test structures
    results = []
    all_plifs = {}  # {species_name: merged_plif_df} — used for heatmap

    for pdb in os.listdir(plip_results_dir):
        if pdb.endswith("_protonated.pdb") and pdb != ref_pdb:
            pdb_name = pdb.replace("_protonated.pdb", "")
            xml = pdb_name + ".xml"
            
            print(f"\nProcessing {pdb_name}...")
            
            # Parse test PLIP data
            test_plip_df = parse_plip_xml(xml)
            
            # Get test van der Waals contacts
            test_vdw_df = get_vdw_contacts(pdb, ligand_name)
            
            # Combine test PLIP and VDW data
            test_df = pd.concat([test_plip_df, test_vdw_df]).reset_index(drop=True)
            
            print(f"  Test PLIF contains {len(test_df)} interactions")
            print(f"    - PLIP-detected: {len(test_plip_df)}")
            print(f"    - VDW contacts: {len(test_vdw_df)}")
            
            # Merge reference and test PLIFs
            merged_plif = merge_plifs(ref_df, test_df)
            all_plifs[pdb_name] = merged_plif

            # Calculate Tanimoto coefficient
            tanimoto = calculate_tanimoto(merged_plif)
            
            shared = merged_plif[(merged_plif["ref"] == 1) & (merged_plif["test"] == 1)].shape[0]
            print(f"  Shared interactions: {shared}")
            print(f"  Tanimoto similarity: {tanimoto}")
            
            # Save individual PLIF comparison
            output_file = f"{pdb_name}_vs_{ref_pdb_name}_PLIF.txt"
            with open(output_file, "w") as f:
                f.write(merged_plif.to_string())
            
            # Add to results
            results.append({
                "test_structure": pdb_name,
                "reference_structure": ref_pdb_name,
                "tanimoto_similarity": tanimoto,
                "ref_interactions": len(ref_df),
                "test_interactions": len(test_df),
                "shared_interactions": shared
            })
    
    # Create results summary table
    print("\n" + "=" * 70)
    print("Creating summary table...")
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("tanimoto_similarity", ascending=False)
    
    # Save summary table — to results/ if provided, otherwise alongside PLIF_files
    summary_dest = summary_dir if summary_dir else plip_results_dir
    os.makedirs(summary_dest, exist_ok=True)
    summary_path = os.path.join(summary_dest, "plif_similarity_summary.csv")
    results_df.to_csv(summary_path, index=False)
    print(f"\nSummary table saved to: {summary_path}")

    # Generate heatmap — one row per species, keeping the best-Tanimoto model
    if all_plifs:
        best_by_species = {}
        for entry in results:
            species = entry["test_structure"].split("_")[0]
            if (species not in best_by_species or
                    entry["tanimoto_similarity"] > best_by_species[species]["tanimoto_similarity"]):
                best_by_species[species] = entry
        best_plifs = {
            e["test_structure"]: all_plifs[e["test_structure"]]
            for e in best_by_species.values()
        }
        tanimoto_scores = {
            e["test_structure"]: e["tanimoto_similarity"]
            for e in best_by_species.values()
        }
        generate_heatmap(best_plifs, ref_df, ref_pdb_name, summary_dest, tanimoto_scores)

    # Organize per-comparison PLIF text files into PLIF_files/
    print("\n" + "=" * 70)
    print("Organizing output files...")
    print("=" * 70)

    output_dir = os.path.join(plip_results_dir, "PLIF_files")
    os.makedirs(output_dir, exist_ok=True)

    files_moved = 0
    for file in os.listdir(plip_results_dir):
        if file.endswith("_PLIF.txt"):   # summary CSV is already saved separately
            try:
                shutil.move(os.path.join(plip_results_dir, file),
                            os.path.join(output_dir, file))
                files_moved += 1
            except shutil.Error:
                print(f"  Warning: {file} already exists in PLIF_files/")

    print(f"Moved {files_moved} PLIF file(s) to PLIF_files/")
    
    return results_df


def main():
    """
    Main integrated workflow for PLIP analysis and PLIF generation.
    
    Workflow:
    1. User provides directory with PDB files (reference + test structures)
    2. PLIP analysis runs on all structures, generating XML reports
    3. Protonated structures are processed with OpenBabel
    4. Results are organized into plip_results directory
    5. User specifies ligand ID and confirms reference structure
    6. PLIFs are generated for all test structures vs reference
    7. Summary table with Tanimoto similarities is created
    8. Individual comparison files are organized for detailed inspection
    
    Benefits of integration:
    - Reduces manual steps and potential for user error
    - Maintains consistent directory structure throughout workflow
    - Provides continuous progress feedback across both stages
    - Ensures all intermediate files are properly organized
    - Allows for error handling across the entire pipeline
    """
    check_tools(["plip", "obabel"])
    print("\n" + "=" * 70)
    print("PLIP ANALYSIS AND PLIF GENERATION PIPELINE")
    print("=" * 70)
    print("\nThis script will:")
    print("1. Run PLIP analysis on all PDB models in your directory")
    print("2. Generate protein-ligand interaction fingerprints (PLIFs)")
    print("3. Compare all test structures against the reference structure")
    print("4. Calculate Tanimoto similarity coefficients")
    print("5. Create organized output files for analysis")

    from utils import resolve_project_dir, load_config, save_config, get_project_paths

    project_dir = resolve_project_dir()
    config = load_config(project_dir)
    paths = get_project_paths(project_dir)

    # ── Models directory ──────────────────────────────────────────────────────
    dir_path = paths["models"]
    if not os.path.isdir(dir_path):
        print(f"\nWarning: models/ not found at {dir_path}")
        dir_path = input("Enter the path to the folder containing docked PDB models: ").strip().strip('"')
        while not os.path.exists(dir_path):
            print("Error: Directory does not exist.")
            dir_path = input("Please enter a valid directory path: ").strip().strip('"')
    else:
        print(f"\nUsing models directory: {dir_path}")

    # Run PLIP analysis — save results into results/PLIP_results/
    os.makedirs(paths["results"], exist_ok=True)
    plip_results_dir = run_plip_analysis(
        dir_path,
        results_dir=os.path.join(paths["results"], "PLIP_results"),
    )

    # ── Ligand residue name ───────────────────────────────────────────────────
    ligand = config.get("ligand_resname")
    if ligand:
        print(f"\nUsing ligand ID from config: {ligand}")
    else:
        ligand = input("\nEnter the ligand ID as it appears in the PDB models (e.g., 'UNL', 'LIG'): ")
    
    # Generate PLIFs — write summary CSV directly to results/
    results_df = generate_plifs(plip_results_dir, ligand,
                                summary_dir=paths["results"])
    
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE")
    print("=" * 70)
    print(f"\nResults are organized in: {plip_results_dir}")
    print(f"  - XML reports and protonated PDBs in main directory")
    print(f"  - PLIF comparisons and summary in PLIF_files subdirectory")


if __name__ == "__main__":
    main()