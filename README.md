# Cross-Species Molecular Docking

## About

Molecular docking is commonly employed to screen large chemical libraries for potential drug candidates. This pipeline repurposes that approach to assess receptor binding across species rather than across chemicals — analogous to "reverse docking," but applied to orthologs of a single receptor type derived from various species. The goal is to generate or support species susceptibility predictions for a chemical of interest by comparing docked poses against those of an empirically derived reference structure.

The final output of this pipeline is a ranked susceptibility assessment in which each test species is assigned a confidence level (Strong, Moderate, or Weak) based on how closely its docked poses resemble those of the reference species across four docking metrics: binding affinity, ligand RMSD, protein-ligand interaction fingerprint (PLIF) Tanimoto similarity, and binding pocket similarity (PPS-score). Ranking is determined by Mahalanobis distance from the reference species centroid in multivariate metric space.

## Citation

If using this repository or any of the code in your work, please cite both:

- Schumann et al., 2024 — original framework: https://doi.org/10.1016/j.comtox.2024.100319
- Schumann et al., 2026 — describes the updated docking methodology: (manuscript in progress)

## Usage

Each script is designed to be run directly from a command-line interface. The simplest approach is to navigate to the `scripts` directory in the cloned repository and call scripts by name:

```
cd path\to\cross-species-docking\scripts
python script_name.py
```

Each script will prompt the user for any required inputs at runtime.

---

## Requirements

### Operating System

This code was developed and tested on Windows 11. Compatibility with other operating systems is not guaranteed.

### 1. Python 3

Download and install [Python 3](https://www.python.org/downloads/).

### 2. Anaconda or Miniconda

[Miniconda](https://docs.anaconda.com/miniconda/install/) is recommended. All scripts must be run within the conda environment described in the Getting Started section below.

### 3. PPS-align

Download the PPS-align source code from the [Zhang Group](https://zhanggroup.org/PPS-align/download.html) and compile it. Compilation requires a g++ compiler-driver; on Windows, [Mingw-w64](https://www.mingw-w64.org/downloads/) or [MSYS2](https://www.msys2.org/) can be used if one is not already installed.

After compiling, add the `PPSalign` executable to your system PATH. Instructions for setting a PATH variable on Windows can be found [here](https://www.java.com/en/download/help/path.html).

### 4. AutoDock Vina

Download the latest release from the [AutoDock Vina GitHub](https://github.com/ccsb-scripps/AutoDock-Vina/releases). Rename the downloaded executables as follows, then add their containing directory to your PATH:

- `vina_1.2.#_win.exe` → `vina.exe`
- `vina_split_1.2.#_win.exe` → `vina_split.exe`

### 5. MUSCLE

Download the [MUSCLE](https://github.com/rcedgar/muscle/releases/tag/5.1.0) executable, rename it `muscle.exe`, and add it to your PATH.

### 6. Meeko

Due to dependency conflicts between Meeko and Python 3.12, installation from source is strongly recommended. After activating the conda environment (see Getting Started), run:

```
git clone https://github.com/forlilab/Meeko.git
cd Meeko
git checkout develop
pip install .
pip install scipy rdkit gemmi tqdm
```

### 7. TM-align

Download and compile TM-align following the instructions on the [TM-align page](https://aideepmed.com/TM-align/). Use a Mingw-w64 or MSYS2 terminal on Windows. Add the compiled executable to your PATH.

### 8. PyMOL (optional)

PyMOL is included in the conda environment and will be available for use within that environment after setup. A standalone installation is only necessary if you want to access PyMOL outside of the `cross-species-docking` environment. Instructions for installing open-source PyMOL on Windows can be found [here](https://pymolwiki.org/index.php/Windows_Install#Open-Source_PyMOL).

---

## Getting Started

**1.** Open an Anaconda Prompt (on Windows, search for "Anaconda" in the taskbar).

**2.** Navigate to the directory where you would like to clone the repository:
```
cd C:\Users\your_username\Documents\docking_analysis
```

**3.** Clone the repository and navigate into it:
```
git clone https://github.com/pschumann4/cross_species_docking.git
cd cross_species_docking
```

**4.** Install `conda-lock`, then create and activate the environment:
```
conda install -c conda-forge conda-lock
conda-lock install --name cross-species-docking conda-lock.yml
conda activate cross-species-docking
```

To use these scripts in the future, activate the environment first with:
```
conda activate cross-species-docking
```

**Note:** It is recommended to create a separate working directory for each analysis. This directory should contain all protein structures and docking input/output files.

---

## Pipeline Overview

The pipeline proceeds in 12 steps, from protein structure acquisition through final susceptibility assessment. It is recommended to update your working directory to `~\cross-species-docking\scripts` before beginning.

**Project configuration:** Each script automatically locates your working directory by searching for a `project_config.json` file. The first time a script runs, you will be prompted to confirm or enter the project directory path; this value is saved to `project_config.json` and reused by all subsequent scripts. Other settings (ligand name, reference PDB path, docking parameters) are saved to this file in the same way, so most prompts only appear once.

---

### Step 1: Obtain protein structures

**Reference structure**

The chemical to be evaluated must have an empirically solved co-crystal structure with the target receptor. For example, butylparaben was crystallized in complex with ESR1 ([PDB: 4MG9](https://www.rcsb.org/structure/4MG9)). Use `smiles_to_pdb.py` to identify candidate PDB entries for your chemical of interest:

```
python smiles_to_pdb.py
```

This outputs a CSV of matching PDB entries with metadata including resolution. Select the highest-resolution structure that contains your ligand of interest bound to the receptor. Avoid using structures with experimental mutations in the binding pocket.

**Test structure set**

A set of species-specific receptor structures is required for cross-species comparison. For a stronger weight of evidence, the [SeqAPASS tool](https://seqapass.epa.gov/seqapass/) (Level 4, I-TASSER homology modeling) is recommended.

Alternatively, use `get_test_strucs.py` to automatically download AlphaFold structural predictions for a diverse set of species. You will need the UniProt ID for your reference protein (e.g., `P10275` for human androgen receptor). The script queries [OrthoDB](https://www.orthodb.org/) to identify orthologs, resolves their UniProt IDs, and downloads available AlphaFold structures along with associated metadata:

```
python get_test_strucs.py
```

This produces three outputs: `orthodb_orthologs.csv`, `alphafold_metadata.csv`, and an `af_structures/` directory containing the downloaded PDB files.

After downloading, manually inspect the structure set for cases where a single species has multiple predicted structures. In these cases, retain the structure with the highest pLDDT confidence score and most relevant description.

**REQUIRED — File naming convention**

PDB file names must follow this format exactly before proceeding:
```
SpeciesName_ProteinSymbol.pdb
```
Examples:
```
Human_AR.pdb
Homosapiens_AR.pdb
Homo-sapiens_AR.pdb
```
Downstream scripts parse species names and protein symbols from these filenames. Files that do not conform to this convention will cause errors in later steps. Rename all files before continuing.

---

### Step 2: Filter structures by structural similarity

Add the reference PDB file to the folder containing the downloaded AlphaFold structures. Then run `run_tmalign.py` to compute TM-scores for all test structures relative to the reference:

```
python run_tmalign.py
```

This outputs `tm_scores.csv` and `tmscore_density.png`. Structures with a TM-score normalized by reference length (TM-score_norm_ref) below 0.5 are considered to have poor structural similarity and should be removed before proceeding. TM-score interpretation:

| TM-score | Interpretation |
|----------|---------------|
| < 0.17 | Random structural similarity |
| ~0.5 | Likely same fold |
| > 0.5 | Likely same topology |
| 1.0 | Identical structures |

---

### Step 3: Align protein structures

Ensure the reference structure is in the same directory as the test structures.

Run `multiple_prot_align.py` to perform a multiple sequence alignment (via MUSCLE) and update residue numbering across all structures, followed by structural alignment via PyMOL. This script also removes extraneous chains and duplicate atoms, and trims each structure to the region aligned with the reference (plus a 10-residue buffer at each terminus):

```
python multiple_prot_align.py
```

Outputs include a set of modified PDB files, a `results/` folder with alignment metadata, and `residue_positions.csv` mapping new residue positions back to their original numbering.

---

### Step 4: Prepare structures with molecular dynamics simulation

Run `mds_structure_prep.py` to clean each structure using PDBFixer (adding missing atoms and resolving common structural issues), then run a short molecular dynamics simulation using OpenMM with the AMBER14 force field to equilibrate each structure and identify flexible binding pocket residues:

```
python mds_structure_prep.py
```

Equilibration is detected automatically using changepoint analysis (PELT algorithm). The recommended simulation length is 100–500 ps; shorter simulations run faster but may not fully capture conformational flexibility. You will also have the option to extract an ensemble of structures from the equilibrated plateau region, which can improve docking accuracy by sampling a range of conformational states.

Outputs include fixed PDB files, per-structure RMSF plots with flexible residues highlighted, a summary of flexible residues for each structure, and (if ensemble extraction is enabled) a set of equilibrated structures.

---

### Step 5: Prepare receptor files for docking

Run `prep_pdbqt.py` to parameterize each receptor and generate rigid and flexible PDBQT files for use with AutoDock Vina:

```
python prep_pdbqt.py
```

---

### Step 6: Generate AutoDock Vina configuration files

Run `generate_config_files.py` to create a Vina configuration file for each receptor:

```
python generate_config_files.py
```

The script auto-resolves the `pdbqt_files/` directory from the project config. Gridbox dimensions are extracted automatically from the reference PDB file and ligand residue name (both read from `project_config.json`); a manual fallback is available if needed. Flexible residue assignments are loaded automatically from `results/flex_residues.txt`. Docking parameters are read from the project config if previously entered, or prompted once and saved.

Recommended settings:

| Parameter | Recommended Value |
|-----------|------------------|
| Number of modes | 5 |
| Energy range | 3 |
| Exhaustiveness | 16 |

---

### Step 7: Run docking simulations and generate models

Ensure the prepared ligand PDBQT file is in the `pdbqt_files/` folder, then run:

```
python run_vina_batch.py
```

This script handles the complete post-docking workflow in one step: it runs AutoDock Vina for all configuration files, selects the pose with the lowest RMSD relative to the reference ligand (using the Hungarian algorithm for optimal atom matching), converts the selected pose to PDB format, and writes a summary of results.

Outputs:
- `docking_results/` — Vina output PDBQT files for the selected pose per receptor
- `docking_results/models/` — combined PDB models (protein + ligand) for each docked complex
- `results/docking_scores.csv` — binding affinities and ligand RMSD values for all receptors

---

### Step 8: Calculate binding pocket similarity scores

Run `run_ppsalign.py` to extract binding pocket structures from each docked complex and compute PPS-scores (binding pocket similarity scores) relative to the reference using PPSalign:

```
python run_ppsalign.py
```

The models directory is resolved automatically from the project config. Output is written to `results/PPS_files/`.

---

### Step 9: Calculate protein-ligand interaction fingerprints (PLIFs)

Run `plip_plif.py` to identify protein-ligand interactions and compute PLIF Tanimoto similarity scores. This uses the [PLIP tool](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index) to detect hydrogen bonds, hydrophobic contacts, and other interaction types, and a custom distance-based algorithm to detect van der Waals interactions. A PLIF is generated for each docked model, and Tanimoto similarity to the reference structure PLIF is calculated:

```
python plip_plif.py
```

The models directory is resolved from the project config. The PLIF Tanimoto summary (`plif_similarity_summary.csv`) is written to `results/`.

---

### Step 10: Generate summary report

Run `get_summary.py` to consolidate all docking metrics into a single summary file:

```
python get_summary.py
```

No file path prompts are issued. All inputs are resolved automatically from the `results/` directory:
- `results/docking_scores.csv` — binding affinities and ligand RMSD (required; generated by Step 7)
- `results/PPS_files/` — PPS-score data (optional; generated by Step 9)
- `results/plif_similarity_summary.csv` — PLIF Tanimoto values (optional; generated by Step 10)

Output: `results/<ligand>_summary.csv`. This file is the sole input to the final susceptibility analysis step. Because RMSD is read directly from `docking_scores.csv`, Step 8 does not need to be run before this step.

---

### Step 11: Perform susceptibility analysis

Run `susceptibility_analysis.py` using the summary report generated in the previous step:

```
python susceptibility_analysis.py
```

The summary CSV is located automatically from `results/`. You will be prompted to specify the name of your reference species (or PDB ID, if the filename was not updated) and whether to apply permissive thresholding (see below).

Output is written to `results/susceptibility_analysis/` and includes:
- `species_summary.csv` — per-species confidence levels and Mahalanobis distances
- `per_model_results.csv` — per-model metric values and pass/fail status
- `pca_plot.png` — PCA visualization of metric space
- `outliers_removed.csv` — models excluded from analysis as outliers
- `analysis_summary.txt` — plain-text summary of results

**Confidence thresholds**

Each species is evaluated against the following default threshold values:

| Metric | Default Threshold |
|--------|-------------------|
| Binding affinity (kcal mol⁻¹) | ≤ −6.0 |
| Ligand RMSD (Å) | ≤ 2.0 |
| PLIF Tanimoto similarity | ≥ 0.5 |
| Protein–pocket similarity score (PPS-score) | ≥ 0.5 |

Confidence levels are assigned based on how many metrics meet or exceed their thresholds:

| Confidence Level | Metrics Passing | Interpretation |
|------------------|-----------------|----------------|
| Strong | 3 – 4 | Docking pose is consistent with the reference across the majority of evaluated criteria; susceptibility is supported with high confidence |
| Moderate | 2 | Partial agreement with the reference; susceptibility is plausible but not fully supported |
| Weak | 0 – 1 | Docking pose diverges substantially from the reference; susceptibility is not well supported |

Species are ranked by Mahalanobis distance from the reference species centroid in multivariate metric space, providing an overall measure of docking similarity that accounts for correlations among metrics.

**Permissive thresholding**

An optional permissive thresholding mode is available, in which the default threshold values are relaxed based on the observed variability within the reference species ensemble. This approach is intended for situations where reference variability is considered high enough that the default thresholds may produce false negatives — that is, where failing a threshold may reflect reference-level noise rather than a genuine difference in binding. It's recommended that permissive thresholding is applied only when there is justification for doing so, as it increases the risk of false positives.

---