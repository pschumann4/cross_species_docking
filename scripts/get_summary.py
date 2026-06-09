"""
Assemble the final per-model summary CSV from all pipeline outputs.

All data sources are resolved automatically from the project's results/ directory:
    results/docking_scores.csv          — binding affinities + lig_rmsd  (required)
    results/PPS_files/*_PPS.txt         — PPS-Score values                (optional)
    results/plif_similarity_summary.csv — PLIF Tanimoto                   (optional)

lig_rmsd is read directly from the selected_mode_rmsd_A column of docking_scores.csv.
run_vina_batch.py already computes this value (Hungarian-algorithm RMSD vs reference)
while selecting the best pose, so running ligand_rmsd.py first is not required.
ligand_rmsd.py remains useful only if you want the RMSD histogram (ligand_rmsd.png).

Output: results/<ligand>_summary.csv
No user prompts are issued — missing optional sources are skipped with a warning.
"""

import os
import sys
import pandas as pd


def get_summary():
    """
    Build the summary CSV from files already present in results/.
    Missing optional sources are reported and left as NaN in the output.
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from utils import resolve_project_dir, load_config, save_config, get_project_paths

    project_dir = resolve_project_dir()
    config      = load_config(project_dir)
    paths       = get_project_paths(project_dir)
    details     = paths["results"]

    print("\n" + "=" * 60)
    print("ASSEMBLING PIPELINE SUMMARY")
    print("=" * 60)

    # ── Ligand name ───────────────────────────────────────────────────────────
    ligand = config.get("ligand_name")
    if not ligand:
        ligand = input("Ligand name not found in config. Enter ligand name: ").strip()
        config["ligand_name"] = ligand
        save_config(project_dir, config)
    print(f"Ligand: {ligand}")

    # ── 1. Docking scores (required) ─────────────────────────────────────────
    scores_csv = os.path.join(details, "docking_scores.csv")
    if not os.path.exists(scores_csv):
        print(f"\nERROR: docking_scores.csv not found in:\n  {details}")
        print("Re-run run_vina_batch.py to regenerate it.")
        return

    scores_df = pd.read_csv(scores_csv).sort_values("protein").reset_index(drop=True)
    print(f"\n[1/4] Docking scores loaded — {len(scores_df)} protein(s)")

    # Parse species and ensemble number from protein name.
    # Examples:
    #   "Chicken_AR_modified_ensemble_1" → species="Chicken", ensemble=1
    #   "Human_AR_modified"              → species="Human",   ensemble=0
    species   = [row.split("_")[0] for row in scores_df["protein"]]
    ensembles = [
        int(row.split("_")[-1]) if row.split("_")[-1].isdigit() else 0
        for row in scores_df["protein"]
    ]

    # lig_rmsd comes directly from docking_scores.csv (selected_mode_rmsd_A).
    # run_vina_batch.py already computes this via the Hungarian algorithm while
    # selecting the best-RMSD pose; ligand_rmsd.py would recompute the same
    # number from the saved model PDB, so there is no reason to run it first.
    summary_df = pd.DataFrame({
        "binding_model":    range(1, len(scores_df) + 1),
        "species":          species,
        "ensemble":         ensembles,
        "binding_affinity": scores_df["best_affinity_kcal_mol"].tolist(),
        "lig_rmsd":         scores_df["selected_mode_rmsd_A"].tolist(),
        "ppsscore":         None,
        "plif_tanimoto":    None,
    })
    print(f"      lig_rmsd pulled from selected_mode_rmsd_A column")

    # ── 2. PPS-Scores (optional) ──────────────────────────────────────────────
    pps_dir = os.path.join(details, "PPS_files")
    if os.path.isdir(pps_dir):
        pps_files = sorted(f for f in os.listdir(pps_dir) if f.endswith("_PPS.txt"))
        ppsscores = []
        for fname in pps_files:
            with open(os.path.join(pps_dir, fname)) as f:
                lines = f.readlines()
            try:
                # PPS-Score is on the third line (index 2), third whitespace-separated token
                ppsscores.append(lines[2].split()[2])
            except (IndexError, ValueError) as e:
                print(f"  Warning: could not parse PPS score from {fname}: {e}")
                ppsscores.append(None)
        summary_df["ppsscore"] = ppsscores
        print(f"[2/3] PPS-Scores loaded — {len(ppsscores)} file(s)")
    else:
        print(f"[2/3] PPS_files/ not found in results/ — ppsscore column left blank")

    # ── 3. PLIF Tanimoto (optional) ───────────────────────────────────────────
    plif_file = os.path.join(details, "plif_similarity_summary.csv")
    if os.path.exists(plif_file):
        plif_df = (
            pd.read_csv(plif_file)
              .sort_values("test_structure")
              .reset_index(drop=True)
        )
        summary_df["plif_tanimoto"] = plif_df["tanimoto_similarity"].values
        print(f"[3/3] PLIF Tanimoto loaded — {len(plif_df)} value(s)")
    else:
        print(f"[3/3] plif_similarity_summary.csv not found in results/ — plif_tanimoto column left blank")

    # ── Write output ──────────────────────────────────────────────────────────
    os.makedirs(details, exist_ok=True)
    output_file = os.path.join(details, f"{ligand}_summary.csv")
    summary_df.to_csv(output_file, index=False)

    print(f"\n{'=' * 60}")
    print(f"Summary saved → {output_file}")
    print("=" * 60)
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    get_summary()
