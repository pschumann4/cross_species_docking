# Cross-Species Molecular Docking
## About
Molecular docking is commonly used to screen large lists of chemicals as potential drug candidates. In this case, molecular docking is employed to screen receptors for making or supporting predictions of species susceptibiilty to chemical effects. This is similar to a "reverse-docking" approach, only instead of testing binding to various receptors, binding is tested against different versions of one receptor type derived from various species (i.e., orthologs). The functions provided in this repository are designed to generate or support existing species susceptibility calls to chemical effects through comparisons with empirically derived reference structures in complex with chemicals of interest.

## Citation
Please cite the use of the code in this repository using the associated article: https://doi.org/10.1016/j.comtox.2024.100319

## Usage
Each script is intended to be ran directly in a command-line interface. The easiest implementation is to simply drag and drop the script into a terminal window after typing "python" or "py".
### Basic command
```
python "path/to/script_file.py"
```

## Requirements

### Operating System
This code was developed using Windows 11. It is unknown to the authors whether these scripts will work as intended on other operating systems.

### 1. Python 3
The user will need to have [Python 3](https://www.python.org/downloads/) installed on their computer.

### 2. Anaconda or Miniconda
[Miniconda](https://docs.anaconda.com/miniconda/install/) is the suggested Anaconda Distribution, but users can use Anaconda as well. All the scripts in this repository will need to be run in your conda environment.

## 3. PPS-align

The PPS-align source code will need to be [downloaded](https://zhanggroup.org/PPS-align/download.html) and compiled.

**NOTE:** To complie PPS-align, you will need a g++ compiler-driver installed. For Windows, you can use [Mingw-w64](https://www.mingw-w64.org/downloads/) or [MSYS2](https://www.msys2.org/), if you don't have one already.

Once installed, "PPSalign" will then need to be set as an environment variable in PATH to run the "run_ppsalign.py" script. Here are instructions on how to set a PATH system variable: https://www.java.com/en/download/help/path.html

### 4. AutoDock Vina
The latest release of AutoDock Vina can be downloaded from [here](https://github.com/ccsb-scripps/AutoDock-Vina/releases).

The "vina_1.2.#_win.exe" and "vina_split_1.2.#_win.exe" will need to be updated to "vina.exe" and "vina_split.exe.", respectively. The directory housing these execuetables will then need to be added to your PATH.

### 5. MUSCLE
The [MUSCLE](https://drive5.com/muscle5/) (MUltiple Sequence Comparison by Log- Expectation) execuetable program will need to be downloaded and added to PATH renamed as "muscle.exe." The program can be downloaded [here](https://github.com/rcedgar/muscle/releases/tag/5.1.0).

### 6. PyMOL (open-source)
Instructions for installing open-source PyMOL for Windows can be found [here](https://pymolwiki.org/index.php/Windows_Install#Open-Source_PyMOL).

## Getting Started
1. Open an Anaconda Prompt (on Windows, go to your search bar and type "Anaconda" and if installed properly, you should see an option to open an Anaconda Prompt).

2. Navigate to a directory where you would like to have the repository cloned. For example:
```
cd C:\Users\pschuman\Documents
```

3. Clone the repository:
```
git clone https://github.com/pschumann4/cross_species_docking.git
```

4. Set up conda environment:
```
conda env create -f cross-species-docking.yml
conda activate cross-species-docking.yml
```

**NOTE:** It is recommended to create a separate folder/directory to run this analysis within. This folder is where you should store all the protein structures and docking input/output files. 

## Steps for cross-species docking analysis

**IMPORTANT** You will need to perform this docking analysis twice -- on the ensemble set and the test set. It is highly recommended that you perform the analysis on the ensemble set *first* and then perform it on the test set. Also note that not all steps apply to each set. Please follow the instructions carefully.

### 1. Get protein structures
- **Reference protein structure**
    The chemical you would like to evaluate for species susceptibility will need to have an empirically solved protein complex. For example, butylparaben was crystallized in complex with ESR1 [PDB: 4MG9](https://www.rcsb.org/structure/4MG9). Download the relevant protein complex for your analysis in PDB format.

- **Ensemble protein set**
    You will need a set of empirically derived structures for the protein you are evaluating. For example, if your chemical is bound to ESR1, you will need to download a set of (ideally unmutated) ESR1 structures from the same species. You can obtain these via an "Advanced Search" within the [RCSB PDB](https://www.rcsb.org/).
    **NOTE**: If downloading a batch of structures from the RCSB PDB, the files might be in ".ent" format, in which case the "multiple_prot_align.py" script will automatically convert these into PDBs. All other formats will need to be converted to PDB prior to performing this analysis.

- **Test protein set**
    A set of species-specific protein structures will need to be assembled. Any combination of sources could be used (e.g., AlphaFold, RCSB PDB, I-TASSER, etc.), although we recommend using the [Sequence Alignment to Predict Across Species Susceptibility (SeqAPASS) tool](https://seqapass.epa.gov/seqapass/) via Level 4 I-TASSER homology modeling for building susceptibility predictions.

    **IMPORTANT**: The PDB file _names_ must be formatted as follows:
    ```
    species_proteinsymbol.pdb
    ```
    Examples: 
    ```
    Human_AR.pdb
    ```
    or
    ```
    Homosapiens_AR.pdb
    ```
    or
    ```
    Homo-sapiens_AR.pdb
    ```

All protein structures must be in PDB format.

### 2. Get chemical structures
You will need to obtain 3D chemical structures from sources like [PubChem](https://pubchem.ncbi.nlm.nih.gov/) in .sdf format. Once downloaded, prepare the chemicals using Meeko by running:
```
mk_prepare_ligand -i molecule.sdf -o molecule.pdbqt
``` 

### 2. Prepare the protein structures
Ensure that the ensemble structures and the test structures are stored in separate directories and that the reference structure is added to both.

**Ensemble set**
Run the "multiple_prot_align.py" script.

In addition to generating a new set of modified PDBs, this will create a folder called "details" in the PDB file directory with information on the alignments as well as a CSV file called "residue_positions.csv" that can be useful for converting the new residue positions back to their original positions, if desired.

**Test set**
Run the "prep_test_receptors.py" script.

You will have the option of running a short molecular dynamics simulation (MDS) to equilibrate the structure using [OpenMM](https://openmm.org/).
This will add significant computational time if you choose to do so. The script will automatically select a structure from the MDS trajectory to estimate the equilibrated system.

Next, run the "multiple_prot_align.py" script on the prepared structures.

### 3. Determine the grid box area for the docking simulation

**NOTE:** Steps 4 - 8 apply to both the ensemble set and test set analyses.

Run the "get_pocket.py" script and when prompted, input "single" and provide the path to the directory containing the aligned reference structure, which should now have the prefix "ref_".
For example, if your reference structure was 4MG9, the aligned reference should now be called "ref_4mg9_modified.pdb".

When asked "Output the gridbox coordinates?" choose "y". This will create a .txt file in the folder called "details" that contains the gridbox information needed for docking.

### 4. Determine flexible residues
Run the "get_flex_residues.py" script using the folder containing the modified (i.e., aligned) structures as the input when prompted.

A .txt file listing the flexible residues for each structure will be saved to the "details" folder.

### 5. Prep structures for docking
Run the "prep_pdbqt.py" script to parameterize each receptor file and generate rigid and flexible PDBQT files.

### 6. Create AutoDock Vina configuration files
Run the "get_config_files.py" script using the "pdbqt_files" folder as an input when prompted.

You will also need to provide the file path information for the gridbox and flexible residues text files, which can be added by dragging and dropping the file into the terminal window at the appropriate prompt.

Recommended configurations:
```
Spacing = 1
Number of modes = 5
Energy range = 10
Exhaustiveness = 16
```

This will generate a configuration file for each receptor.

### 7. Perform the docking simulation
Add the chemical PDBQT file to the "pdbqt_files" folder.

Then, run the "run_vina_batch.py" script.

The results will be added to a new folder called "vina_output".

### 8. Generate binding models from Vina outputs
Copy and paste the aligned receptor PDB files (ending with "_modified.pdb") into the "vina_output" folder including the aligned reference structure.

Then, run the "generate_models.py" script.

This will combine all corresponding flexible residues, rigid residues, and ligand poses into a PDB file. These PDB models will be saved to a new folder called "models".

### 9. Calculate ligand RMSD
Run the "ligand_rmsd.py" script using the "models" directory path as an input when prompted.

You will also be asked, "Would you like to filter the best/worst models based on RMSD? (y/n)". If you are performing the ensemble set analysis, then we suggest that you say "y". Sorting the results in this way simplifies the categorization step when generating susceptibility predictions, which are binary -- either "yes" or "no". When sorting the models in this way, you will also be asked to sort the "vina_output", which you should also input "y". Doing this will also reduce computational time on the next steps of the analysis.

### 10. Calculate binding pocket similarity scores
Run the "get_pocket.py" script, again, but this time specifying that you want to run it over a "directory" when prompted. Use the path to the "models" folder as the input (or "filtered_models" if performing the ensemble set analysis).

This will generate a new folder called "binding_sites".

Run the "run_ppsalign.py" script using the "binding_sites" folder as an input.

This will output a folder called "PPS_files" that contains all the PPS Score information.

**NOTE:** With all the folders within folders being generated, it is highly recommended that you organize these output folders into a single location on your computer like the parent directory that you originally created for your receptor files.

### 11. Perform the protein-ligand interaction fingerprint (PLIF) analysis
Run the "run_plip.py" script using the "models" or "filtered_models" (for the ensemble set) directory path.

If PLIP was not installed properly, try running ```conda install -c conda-forge plip```. If issues persist, this analysis can also be performed using the [PLIP web-based tool](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index). Save the .xml output and the protonated versions of the inputted PDB. **NOTE:** If done this way, ensure that the protonated PDBs are fully protonated, including the ligand atoms before moving to the next step.

The protonated PDB and the .xml file for each query structure and the reference structure will need to be in the same "plip_results" folder.

Once the PLIP .xml files and protonated PDBs are generated for each model, run the "get_plif.py" script using the "plip_results" folder as the input when prompted.

### 12. Generate a summary report file
Run the "get_summary.py" script. The user will be prompted for all relevant information to generate the summary report file.

### 13. Perform cluster analysis
Once you have a summary report for the ensemble set analysis and the test set analysis, run the "cluster_analysis_kNN.py" script.

**NOTE:** A self-docking simulation should have been performed automatically using your reference structure. You will need to determine which self-docking result is the "best" before running this script. Typically, this is the model with the lowest calculated ligand RMSD.

The main output of this script is a text file listing the species predicted as "susceptible" or "not likely susceptible" according to these docking analyses.
Briefly, if a query model is found to be within the cluster of the reference structure's self-docking model, then that species would be considered likely susceptible to the effects of the bound chemical.

## Further information
We refer the user to the journal article describing this work in detail for more information: https://doi.org/10.1016/j.comtox.2024.100319

## Disclaimer
The United States Environmental Protection Agency (EPA) GitHub project code is provided on an “as is” basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.
