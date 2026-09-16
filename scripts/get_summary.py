"""
Assemble the final per-model summary CSV from all pipeline outputs.

All data sources are resolved automatically from the project's results/ directory:
    results/docking_scores.csv          — binding affinities + lig_rmsd  (required)
    results/PPS_files/*_PPS.txt         — PPS-Score values                (optional)
    results/plif_similarity_summary.csv — PLIF Tanimoto                   (optional)

binding_affinity and lig_rmsd are read from the binding_affinity_kcal_mol and
lig_rmsd_A columns of docking_scores.csv. Both describe Vina's rank-1 (best-scored)
pose — the same pose the model PDB, PLIF and PPS are derived from — so all four
metrics describe one physical binding event. run_vina_batch.py computes lig_rmsd
(Hungarian-algorithm RMSD vs reference), so running ligand_rmsd.py first is not
required; it remains useful only for the RMSD histogram (ligand_rmsd.png).

Output: results/<ligand>_summary.csv
No user prompts are issued — missing optional sources are skipped with a warning.
"""

import os
import sys
import pandas as pd

# Pipeline-added filename suffixes, stripped (longest first) to recover the
# canonical protein key shared across all three data sources. The protein name
# in docking_scores.csv (e.g. "Chicken_AR_modified_ensemble_1") is the key;
# PLIF rows carry a trailing "_model" and PPS files additionally carry
# "_bindingsite" plus the "_PPS" stem, so those must be removed before joining.
_PROTEIN_KEY_SUFFIXES = ("_bindingsite", "_model")


def canonical_protein_key(name):
    """
    Reduce a PLIF test_structure name or PPS file stem to the canonical protein
    key used in docking_scores.csv by stripping pipeline-added suffixes from the
    end. Idempotent for names that are already canonical.
    """
    key = name
    changed = True
    while changed:
        changed = False
        for suffix in _PROTEIN_KEY_SUFFIXES:
            if key.endswith(suffix):
                key = key[: -len(suffix)]
                changed = True
    return key


def _warn_unmatched(source, protein_keys, value_map):
    """
    Warn about join mismatches between docking_scores proteins and an optional
    source, in both directions:
      - proteins with no value in this source (left as NaN)
      - source entries that matched no protein (silently dropped without this)
    """
    proteins = set(protein_keys)
    missing_for_protein = [k for k in protein_keys if k not in value_map]
    orphan_source_keys  = [k for k in value_map if k not in proteins]
    if missing_for_protein:
        print(f"  ⚠ {source}: {len(missing_for_protein)} protein(s) have no value "
              f"(left blank): {', '.join(missing_for_protein)}")
    if orphan_source_keys:
        print(f"  ⚠ {source}: {len(orphan_source_keys)} entry(ies) matched no protein "
              f"and were ignored: {', '.join(orphan_source_keys)}")


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

    # Require the canonical-pose columns. Older docking_scores.csv files used
    # best_affinity_kcal_mol / selected_mode_rmsd_A, where affinity (mode 1) and
    # RMSD (best-RMSD mode) described different poses; refuse those rather than
    # silently mixing incoherent metrics.
    required_score_cols = {"binding_affinity_kcal_mol", "lig_rmsd_A"}
    missing_score_cols = required_score_cols - set(scores_df.columns)
    if missing_score_cols:
        print(f"\nERROR: docking_scores.csv is missing column(s): {sorted(missing_score_cols)}")
        print("It looks like it was produced by an older run_vina_batch.py.")
        print("Re-run run_vina_batch.py so all four metrics describe the same (rank-1) pose.")
        return

    # Parse species and ensemble number from protein name.
    # Examples:
    #   "Chicken_AR_modified_ensemble_1" → species="Chicken", ensemble=1
    #   "Human_AR_modified"              → species="Human",   ensemble=0
    species   = [row.split("_")[0] for row in scores_df["protein"]]
    ensembles = [
        int(row.split("_")[-1]) if row.split("_")[-1].isdigit() else 0
        for row in scores_df["protein"]
    ]

    # binding_affinity and lig_rmsd both come from the canonical (mode-1) pose,
    # the same pose the model PDB — and thus the PLIF and PPS metrics — are built
    # from, so all four metrics describe one physical binding event.
    #
    # The protein column is the canonical key. PPS and PLIF values are joined
    # onto it BY KEY (not by row position) so that differing sort orders or a
    # protein missing from one optional source can never silently shift a metric
    # onto the wrong species.
    protein_keys = scores_df["protein"].tolist()
    summary_df = pd.DataFrame({
        "binding_model":    range(1, len(scores_df) + 1),
        "species":          species,
        "ensemble":         ensembles,
        "binding_affinity": scores_df["binding_affinity_kcal_mol"].tolist(),
        "lig_rmsd":         scores_df["lig_rmsd_A"].tolist(),
        "ppsscore":         None,
        "plif_tanimoto":    None,
    })
    print(f"      binding_affinity and lig_rmsd taken from the rank-1 pose")

    # ── 2. PPS-Scores (optional) ──────────────────────────────────────────────
    pps_dir = os.path.join(details, "PPS_files")
    if os.path.isdir(pps_dir):
        pps_map = {}
        for fname in sorted(f for f in os.listdir(pps_dir) if f.endswith("_PPS.txt")):
            key = canonical_protein_key(fname[: -len("_PPS.txt")])
            with open(os.path.join(pps_dir, fname)) as f:
                lines = f.readlines()
            try:
                # PPS-Score is on the third line (index 2), third whitespace-separated token
                pps_map[key] = lines[2].split()[2]
            except (IndexError, ValueError) as e:
                print(f"  Warning: could not parse PPS score from {fname}: {e}")
                pps_map[key] = None

        summary_df["ppsscore"] = [pps_map.get(k) for k in protein_keys]
        matched = sum(1 for k in protein_keys if k in pps_map)
        print(f"[2/3] PPS-Scores loaded — {len(pps_map)} file(s), {matched}/{len(protein_keys)} matched to proteins")
        _warn_unmatched("PPS", protein_keys, pps_map)
    else:
        print(f"[2/3] PPS_files/ not found in results/ — ppsscore column left blank")

    # ── 3. PLIF Tanimoto (optional) ───────────────────────────────────────────
    plif_file = os.path.join(details, "plif_similarity_summary.csv")
    if os.path.exists(plif_file):
        plif_df = pd.read_csv(plif_file)
        plif_map = {
            canonical_protein_key(str(ts)): tan
            for ts, tan in zip(plif_df["test_structure"], plif_df["tanimoto_similarity"])
        }
        summary_df["plif_tanimoto"] = [plif_map.get(k) for k in protein_keys]
        matched = sum(1 for k in protein_keys if k in plif_map)
        print(f"[3/3] PLIF Tanimoto loaded — {len(plif_map)} value(s), {matched}/{len(protein_keys)} matched to proteins")
        _warn_unmatched("PLIF", protein_keys, plif_map)
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
