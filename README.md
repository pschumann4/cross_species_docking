# Cross-Species Molecular Docking
## About
Molecular docking is commonly used to screen large lists of chemicals as potential drug candidates. In this case, molecular docking is employed to screen receptors for making or supporting predictions of species susceptibiilty to chemical effects. This is similar to a "reverse-docking" approach, only instead of testing binding to various receptors, binding is tested against one receptor type. Each version of the receptor is derived for different species protein information (i.e., protein orthologs). The functions provided in this repository are designed to derive species susceptibility calls to chemical effects namely through comparisons with empirically derived reference structures.

## Citation
Citation will become available upon journal publication.

## Usage
Each script is intended to be ran directly in a terminal or command line.
### Basic command
```
python "path/to/script_file.py"
```

## Requirements

### Operating System
This code was developed using Windows 11. Therefore, it is unknown to the authors whether these scripts will work as intended in other systems.

### 1. Python 3
The user will need to have [Python](https://www.python.org/downloads/) installed on their computer. 

See requirements.txt for required Python packages.

These can be installed using pip by running the following in the terminal:
```
pip install -r requirements.txt 
```

### 2. MGLTools
MGLTools is used for receptor and ligand prep and will need to be [installed](https://ccsb.scripps.edu/mgltools/downloads/).

MGLTools will need to be installed so that the user can utilize the MolKit module. The user will also need to ensure that Python 2 can be accessed. Depending on the user's OS, the system "python" might be replaced with Python 2 by default after installing MGLTools. 

**NOTE**: If errors are encountered with running the "prep_receptors.py" script, check that "python" is calling v2 in the PATH. Otherwise, update the script directly by replacing the command "python" with whatever variable your system is using to call Python 2.

### 3. PPS-align

The PPS-align source code will need to be [downloaded](https://zhanggroup.org/PPS-align/download.html) and compiled.

**NOTE:** To complie PPS-align, you will need the g++ compiler-driver installed. For Windows, you can use [MSYS2](https://www.msys2.org/), if you don't have one already.

"PPSalign" will then need to be set as an environment variable in PATH to run the "ppsalign_loop.py" script. Otherwise, you can update the script directly by replacing the "PPSalign" variable with the path to the compiled PPS-align program.

### 4. Protein-Ligand Interaction Profiler (PLIP) Tool
It is recommended that this is [installed](https://github.com/pharmai/plip) using conda:
```
conda install -c conda-forge plip
```

The "batch_plipcmd.bat" function should then be ran using a terminal within the Conda environment containing this install.

### 5. AutoDock Vina
The latest release of AutoDock Vina can be installed from [here](https://github.com/ccsb-scripps/AutoDock-Vina/releases).

The variable, "vina", will need to be added to the users PATH to run the "batch_vina.bat" function directly. Otherwise, the user will need to replace the instance of "vina" in this script with the path to the vina.exe found in their system.

### 6. MUSCLE
The [MUSCLE](https://drive5.com/muscle5/) (MUltiple Sequence Comparison by Log- Expectation) execuetable program will need to be downloaded and added to PATH renamed as "muscle.exe." The program can be downloaded [here](https://github.com/rcedgar/muscle/releases/tag/5.1.0).
## Steps for cross-species docking analysis
**IMPORTANT**: To make susceptibility predictions on a set of species structures, the user will need to perform this process twice. Once with an ensemble set of reference species structures, and then again with the set of query species structures.

### 1. Get protein structures
For an ensemble analysis, we recommend the user obtains these via an "Advanced Search" within the [RCSB PDB](https://www.rcsb.org/).

For a query structure analysis, we recommend that the user obtains species protein structures from the [Sequence Alignment to Predict Across Species Susceptibility (SeqAPASS) tool](https://seqapass.epa.gov/seqapass/) via Level 4 evaluations. However, structures could be derived from any available source (i.e., [RCSB PDB](https://www.rcsb.org/) or [AlphaFold DB](https://alphafold.com/)). 

All structure files must be in PDB format! 

**NOTE**: If downloading a batch of structures from the RCSB PDB, the files might be in ".ent" format, in which case the "multiple_prot_align.py" script will automatically convert these into PDBs. All other formats will need to be converted to PDB prior to performing this analysis. 

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
### 2. Multiple protein sequence and structural alignments
Run the "multiple_prot_align.py" script.

In addition to generating a new set of modified PDBs, this will create a folder called "details" in the PDB file directory with information on the alignments as well as a CSV file called "residue_positions.csv" that can be useful for converting the new residue positions back to their original positions, if desired.
### 3. Determine the grid box area for the docking simulation
Run the "get_poc.py" script and when prompted, input "single" and select the reference structure with the chemical of interest bound, which can be found in the "original_structures" folder.

The user will be prompted as to whether they would like to generate a gridbox file. Select "y". We recommend that you save the outputted files in the "details" folder.
### 4. Determine flexible residues
Run the "get_flex_residues.py" script.

We recommend that the outputted "flex_residues.txt" file is saved in the "details" folder.
### 5. Prep structures for docking
Run the "prep_receptors.py" script.

Ensure that "python" is calling Python 2. Otherwise, the user will need to manually edit the scripts to match their system's python PATH variable.

This step is used to generate PDBQT files and any flex/rigid files for a flexible-receptor docking approach.

Occasionally, files will fail to prep and will be moved into a separate folder called "failed_preps" to indicate so. These failed preps will likely need to be preppared using the AutoDock Tools GUI.
### 6. Create AutoDock Vina configuration files
Run the "get_config_files.py" script.

Recommended configurations:
Spacing = 1
Number of modes = 3-5
Energy range = 10
Exhaustiveness = 8

The user will need the grid box information generated from step 3 and the "flex_residues.txt" file generated in step 4.
### 7. Perform the docking simulation
Ensure that the chemical of interest to be used in the docking simulation has been prepped as a PDBQT file with polar hydrogens and charges added. Then add the ligand PDBQT file to the same directory as the protein receptor PDBQT files.

Using AutoDock Vina, the script "batch_vina.bat" can be ran on Windows OS.'

After the docking simulation is complete, run the "batch_vina_split.bat" function to generate separate PDBQT files for each predicted binding mode. This function will move the output files into a new folder called "vina_output."
### 8. Generate binding models from Vina outputs
The PDB files generated from the multiple protein alignment (i.e., the PDB files ending with "_modified.pdb") as well as the reference PDB (i.e., the PDB file starting with "ref") will need to be copied into the "vina_output" folder.

Then, run the "add_ligand.py" script using that folder to add the Vina-generated poses to the protein structures.

This script will move all of the PDB models into a new folder called "models."
### 9. Calculate ligand RMSD
Run the "ligand_rmsd.py" script.

**NOTE**: If the user is performing this step for the ensemble docking phase, we recommend that they choose "y" when prompted to filter the models for the "best" and "worst."
### 10. Calculate binding pocket similarity scores
Run the "get_poc.py" script, but this time specifying "directory" when prompted.

If the user preferred to combine all of the pocket files into one, then they can skip the use of the "ppsalign_loop.py" script and instead run: 
```
PPSalign -is1vs1 N QUERYs.poc TEMPL.poc
```
(replacing "QUERYs.poc" with the combined .poc file and "TEMPL.poc" with the reference .poc file.)

Otherwise, run the "ppsalign_loop.py" script to calculate PPS-scores from the directory of .poc files.
### 11. Perform the protein-ligand interaction fingerprint (PLIF) analysis
Using the environment in which the "plip" module was installed, run the "batch_plipcmd.bat" script. 

If PLIP was not installed or if the user is not working within Windows, this analysis can also be performed using the [PLIP web-based tool](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index). Save the .xml output and the protonated versions of the inputted PDB. **NOTE**: If done this way, ensure that the protonated PDBs are fully protonated, including the ligand atoms before moving to the next step.

Once the PLIP .xml files and protonated PDBs are generated for each model, run the "get_plif.py" script.

The user will need to have the protonated PDB and the .xml file for each query structure and the reference structure in the same folder.
### 12. Generate a summary report file
Run the "get_summary.py" script. The user will be prompted for all relevant information to generate the summary report file.
### 13. Perform cluster analysis
Assuming that the user has a summary report for the ensemble analysis and the query structure analysis, run the "cluster_analysis_kNN.py" script. 
## Interpretation of results
We refer the user to the corresponding journal article (in process) for an in-depth explanation of the results from this analysis.

Briefly, if a query model is found to be within the cluster of the reference structure's self-docking model, then that species would be considered likely susceptible to the effects of the bound chemical.
## Disclaimer
The United States Environmental Protection Agency (EPA) GitHub project code is provided on an “as is” basis and the user assumes responsibility for its use. EPA has relinquished control of the information and no longer has responsibility to protect the integrity , confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by EPA. The EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by EPA or the United States Government.
