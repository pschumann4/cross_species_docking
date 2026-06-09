"""
run_tmalign.py
==============================
TM-align for evaluating structural similarity of test models to a reference structure.
This is useful for filtering out poor models before docking, and for assessing the overall quality of the test set.
The script performs the following steps:
1. Load and preprocess data - identify PDB files in directory
2. Identify reference model (prefixed with "ref_")
3. Run TM-align for each test model against the reference
4. Extract TM-scores (both normalizations) and alignment statistics
5. Output results to CSV
6. Generate KDE density plot with cutoff at 0.5
TM-score interpretation:
- TM-score < 0.17: random similarity
- TM-score ~ 0.5: likely same fold
- TM-score > 0.5: likely same topology
- TM-score = 1: identical structures

Scores normalized by the reference length (TM-score_norm_ref) that are below 0.5 are
generally considered to indicate poor structural similarity, while those above 0.5 suggest a reasonable fold match.

Outputs
-------
- tm_scores.csv: A CSV file containing TM-scores and alignment statistics for each test model
- tmscore_density.png: A KDE density plot of TM-score_norm_ref values
"""

import os
import sys
import subprocess
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import check_tools


def plot_tmscore_density(scores, output_path, cutoff=0.5):
    """
    Generate a KDE density plot of TM-scores with a vertical cutoff line.

    The density is estimated via scipy.stats.gaussian_kde, which uses
    Scott's rule by default to select bandwidth: h = n^(-1/5) * sigma.
    The area under the curve is split at `cutoff` and filled with
    distinct colours to make the proportion above/below immediately legible.

    Parameters
    ----------
    scores      : array-like of float   TM-score_norm_ref values
    output_path : str                   Path for the saved PNG
    cutoff      : float                 Threshold line (default 0.5)
    """
    CUTOFF      = cutoff
    COLOR_FAIL  = "#CA043F"   # red  — TM < 0.5
    COLOR_PASS  = "#0B6FC7"   # blue — TM ≥ 0.5
    COLOR_LINE  = "#FD7D11"   # orange — cutoff
    BG          = "#ffffff"
    TEXT     = "#0a0a0a"
    MUTED       = "#0a0a0a"
    SURFACE        = "#ffffff"

    scores = np.asarray(scores)
    n = len(scores)

    # Fit KDE using Scott's rule bandwidth
    kde = gaussian_kde(scores)
    bw_used = kde.factor * scores.std(ddof=1)

    # Evaluate KDE on a dense grid across [0, 1]
    x_grid = np.linspace(0, 1, 500)
    y_grid = kde(x_grid)

    # Split grid at cutoff for two-colour fill
    mask_fail = x_grid <= CUTOFF
    mask_pass = x_grid >= CUTOFF

    # Proportions: integrate KDE on each side (trapezoidal rule)
    prop_fail = np.trapz(y_grid[mask_fail], x_grid[mask_fail])
    prop_pass = np.trapz(y_grid[mask_pass], x_grid[mask_pass])

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 5), facecolor=BG)
    ax.set_facecolor(SURFACE)

    # Filled areas
    ax.fill_between(x_grid, y_grid, where=mask_fail,
                    color=COLOR_FAIL, alpha=0.55, linewidth=0)
    ax.fill_between(x_grid, y_grid, where=mask_pass,
                    color=COLOR_PASS, alpha=0.55, linewidth=0)

    # KDE curve outline
    ax.plot(x_grid[mask_fail], y_grid[mask_fail], color=COLOR_FAIL, lw=2)
    ax.plot(x_grid[mask_pass], y_grid[mask_pass], color=COLOR_PASS, lw=2)

    # Cutoff vertical line
    ax.axvline(CUTOFF, color=COLOR_LINE, lw=1.5, linestyle="--")
    ax.text(CUTOFF + 0.01, ax.get_ylim()[1] * 0.97,
            f"cutoff = {CUTOFF}", color=COLOR_LINE,
            fontsize=9, va="top", fontfamily="monospace")

    # Proportion annotations inside each region
    ax.text(CUTOFF / 2, max(y_grid) * 0.55,
            f"{prop_fail*100:.1f}%\n(TM < {CUTOFF})",
            color=COLOR_FAIL, fontsize=10, ha="center",
            fontweight="bold", fontfamily="monospace")
    ax.text(CUTOFF + (1 - CUTOFF) / 2, max(y_grid) * 0.55,
            f"{prop_pass*100:.1f}%\n(TM ≥ {CUTOFF})",
            color=COLOR_PASS, fontsize=10, ha="center",
            fontweight="bold", fontfamily="monospace")

    # Labels and formatting
    ax.set_xlabel("TM-score (normalized by reference length)", color=TEXT, fontsize=11)
    ax.set_ylabel("Density", color=TEXT, fontsize=11)
    ax.set_title(
        f"TM-score Distribution  |  n = {n}",
        color=TEXT, fontsize=11, pad=12
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)
    ax.tick_params(colors=MUTED)
    for spine in ax.spines.values():
        spine.set_edgecolor("#30363d")
    ax.grid(False)

    # Legend
    legend = ax.legend(
        handles=[
            mpatches.Patch(color=COLOR_FAIL, alpha=0.7, label=f"TM < {CUTOFF}"),
            mpatches.Patch(color=COLOR_PASS, alpha=0.7, label=f"TM ≥ {CUTOFF}"),
        ],
        facecolor=BG, edgecolor="#30363d", labelcolor=TEXT, fontsize=9
    )

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor=BG)
    plt.close(fig)

    print(f"\nDensity plot saved to: {output_path}")
    print(f"  Above cutoff: {prop_pass*100:.1f}%  ({int(round(prop_pass*len(scores)))} / {len(scores)} models)")
    print(f"  Below cutoff: {prop_fail*100:.1f}%  ({int(round(prop_fail*len(scores)))} / {len(scores)} models)")
    print(f"  KDE bandwidth (Scott's rule): {bw_used:.4f}")
    

def run_tmalign():    
    """
    Calculate TM-scores across model directory using TM-align executable,
    then generate a KDE density plot of TM-score_norm_ref with a 0.5 cutoff line.

    Methodology:
    1. Load and preprocess data - identify PDB files in directory
    2. Identify reference model (prefixed with "ref_")
    3. Run TM-align for each test model against the reference
    4. Extract TM-scores (both normalizations) and alignment statistics
    5. Output results to CSV
    6. Generate KDE density plot with cutoff at 0.5
    """

    check_tools(["TMalign"])
    # ── 1. Directory input ────────────────────────────────────────────────────
    pdb_dir = input("Enter the directory containing the PDB models: ").strip()

    while not os.path.exists(pdb_dir):
        print("The directory does not exist.")
        pdb_dir = input("Enter the directory containing the PDB models: ").strip()

    # ── 2. Identify reference PDB (prefixed "ref_") ───────────────────────────
    ref_pdb = None
    pdb_files = [i for i in os.listdir(pdb_dir) if i.endswith(".pdb")]
    
    for file in pdb_files:
        if file.startswith("ref_"):
            ref_pdb = file
            print(f"Found reference structure: {ref_pdb}")
            break
    
    while ref_pdb is None:
        ref_pdb = input(
            "Please enter the filename of the reference PDB file "
            "(in the specified directory): "
        ).strip()
        if not ref_pdb.endswith(".pdb"):
            ref_pdb += ".pdb"
        if not os.path.exists(os.path.join(pdb_dir, ref_pdb)):
            print(f"The reference PDB file '{ref_pdb}' does not exist in {pdb_dir}.")
            ref_pdb = None

    # ── 3. Run TM-align for each test model against the reference ─────────────
    results = []
    ref_pdb_path = os.path.join(pdb_dir, ref_pdb)
    
    print(f"\nRunning TM-align against reference: {ref_pdb}")
    print("-" * 60)
    
    for file in pdb_files:
        if file == ref_pdb:
            continue

        test_pdb = os.path.join(pdb_dir, file)
        command = f'TMalign "{test_pdb}" "{ref_pdb_path}"'

        try:
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=pdb_dir
            )
            stdout, stderr = process.communicate()

            if process.returncode != 0:
                print(f"Warning: TM-align returned error for {file}")
                print(f"  {stderr.decode()}")
                continue

            # ── 4. Parse TM-align output ──────────────────────────────────────
            output = stdout.decode()

            tm_score_chain1 = None  # normalized by test length
            tm_score_chain2 = None  # normalized by reference length
            rmsd = None
            aligned_length = None
            seq_id = None
            chain1_length = None
            chain2_length = None

            for line in output.splitlines():
                if line.startswith("Length of Chain_1:"):
                    m = re.search(r'Length of Chain_1:\s*(\d+)', line)
                    if m:
                        chain1_length = int(m.group(1))

                elif line.startswith("Length of Chain_2:"):
                    m = re.search(r'Length of Chain_2:\s*(\d+)', line)
                    if m:
                        chain2_length = int(m.group(1))

                elif line.startswith("Aligned length="):
                    m_al  = re.search(r'Aligned length=\s*(\d+)', line)
                    m_rm  = re.search(r'RMSD=\s*([\d.]+)', line)
                    m_sid = re.search(r'Seq_ID=n_identical/n_aligned=\s*([\d.]+)', line)
                    if m_al:  aligned_length = int(m_al.group(1))
                    if m_rm:  rmsd           = float(m_rm.group(1))
                    if m_sid: seq_id         = float(m_sid.group(1))

                elif line.startswith("TM-score=") and "Chain_1" in line:
                    m = re.search(r'TM-score=\s*([\d.]+)', line)
                    if m:
                        tm_score_chain1 = float(m.group(1))

                elif line.startswith("TM-score=") and "Chain_2" in line:
                    m = re.search(r'TM-score=\s*([\d.]+)', line)
                    if m:
                        tm_score_chain2 = float(m.group(1))

            results.append({
                "Test_Model":         file,
                "Chain1_Length":      chain1_length,
                "Chain2_Length":      chain2_length,
                "TM-score_norm_test": tm_score_chain1,
                "TM-score_norm_ref":  tm_score_chain2,
                "RMSD":               rmsd,
                "Aligned_Length":     aligned_length,
                "Seq_Identity":       seq_id,
            })

            if tm_score_chain2 is not None:
                print(f"  {file:50s}  TM-score(ref): {tm_score_chain2:.5f}")
            else:
                print(f"  {file:50s}  [parsing failed]")

        except Exception as e:
            print(f"Error processing {file}: {e}")

    # ── 5. Save CSV ────────────────────────────────────────────────────────────
    if not results:
        print("No results generated. Check input files and TM-align installation.")
        return

    results_df = pd.DataFrame(results)
    output_csv = os.path.join(pdb_dir, "tm_scores.csv")
    results_df.to_csv(output_csv, index=False)

    print("\n" + "=" * 60)
    print(f"TM-align complete. Results saved to: {output_csv}")
    print(f"Total models analysed: {len(results)}")
    print("=" * 60)

    valid_scores = results_df["TM-score_norm_ref"].dropna()
    if len(valid_scores) > 0:
        print(f"\nTM-score summary (normalized by reference):")
        print(f"  Mean:   {valid_scores.mean():.5f}")
        print(f"  Median: {valid_scores.median():.5f}")
        print(f"  Min:    {valid_scores.min():.5f}")
        print(f"  Max:    {valid_scores.max():.5f}")
        print(f"  Std:    {valid_scores.std():.5f}")

    # ── 6. KDE density plot ────────────────────────────────────────────────────
    plot_tmscore_density(
        scores=valid_scores.values,
        output_path=os.path.join(pdb_dir, "tmscore_density.png"),
        cutoff=0.5
    )


if __name__ == "__main__":
    run_tmalign()