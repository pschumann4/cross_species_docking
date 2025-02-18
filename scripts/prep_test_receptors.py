import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
import ruptures as rpt
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Define function that runs the PDBFixer
def run_pdbfixer(pdb_name, output_dir):
    """
    Run PDBFixer on a PDB file to fix common issues
    
    Parameters:
    -----------
    pdb_name : str
        Path to the input PDB file
    output_dir : str
        Directory to save output files
    
    Returns:
    --------
    str: Path to the fixed PDB file
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    print("\nPreparing", os.path.basename(pdb_name), "using PDBFixer...")
    
    # Fix the PDB file
    fixer = PDBFixer(pdb_name)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    
    # Output fixed PDB
    fixed_filename = os.path.join(output_dir, os.path.basename(pdb_name).replace(".pdb", "_fixed.pdb"))
    PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_filename, 'w'))

    return fixed_filename

def run_mds(processed_filename, output_dir, mds_time=10):
    """
    Run Molecular Dynamics Simulation on a PDB file

    Parameters:
    -----------
    processed_filename : str
        Path to the processed PDB file
    output_dir : str
        Directory to save output files
    mds_time : int, optional
        Total simulation time in ns (default: 10)

    Returns:
    --------
    str: Path to the processed PDB file
    str: Path to the MDS output PDB file
    float: Time of equilibration
    str: Path to the equilibrated PDB file
    """
    # Get the pdb_name
    pdb_name = os.path.basename(processed_filename).replace("_fixed.pdb", ".pdb")

    # Perform Molecular Dynamics Simulation
    print("\nPerforming " + str(mds_time) + " ns MD simulation for equilibration...")
    mds_output_name = os.path.join(output_dir, os.path.basename(pdb_name).replace(".pdb", "-mds.pdb"))
    pdb = PDBFile(processed_filename)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield, padding=1.0*nanometer)
    
    system = forcefield.createSystem(modeller.topology, 
                                     nonbondedMethod=PME,
                                     nonbondedCutoff=1*nanometer, 
                                     constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy()
    simulation.reporters.append(PDBReporter(mds_output_name, 1000))
    simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
            potentialEnergy=True, temperature=True, volume=True))
    simulation.step(int(mds_time)*1000)

    # Remove the processed PDB file
    os.remove(processed_filename)

    # Determine approximate equilibration point
    print("\nDetermining approximate equilibration point based on RMSD...\n")
    u = mda.Universe(mds_output_name, dt = 1000.0)
    reference = u.select_atoms("protein")
    R = RMSD(u, reference, select="protein")
    R.run()
    
    time = (R.results.rmsd[:, 1])
    rmsd_values = R.results.rmsd[:, 2]

    # Close the Universe after RMSD analysis
    u.trajectory.close()
    del u
    
    def estimate_plateau_point(rmsd_values, time):
        # Plateau detection using PELT algorithm
        algo = rpt.Pelt(model="rbf").fit(rmsd_values)
        result = algo.predict(pen=1)
        plateau_start_index = result[0]
        
        # Calculate plateau statistics
        plateau_values = rmsd_values[plateau_start_index:]
        plateau_average = np.mean(plateau_values)
        plateau_std = np.std(plateau_values)
        
        # Find the RMSD value closest to the mean
        differences = np.abs(plateau_values - plateau_average)
        closest_to_mean_local_idx = np.argmin(differences)
        closest_to_mean_index = plateau_start_index + closest_to_mean_local_idx
    
        return {
            'start_index': plateau_start_index,
            'start_time': time[plateau_start_index],
            'plateau_average': plateau_average,
            'plateau_std': plateau_std,
            'mean_representative_index': closest_to_mean_index,
            'mean_representative_time': time[closest_to_mean_index]
        }

    # Estimate plateau point
    plateau_results = estimate_plateau_point(rmsd_values, time)
    plateau_index = plateau_results['mean_representative_index']
    plateau_time = plateau_results['start_time']

    # Plot RMSD with enhanced visualization
    plt.figure(figsize=(10,6))
    plt.plot(time, rmsd_values, marker='o', label='RMSD')

    # Add plateau line and indicators
    plt.axvline(x=plateau_time, color='r', linestyle='--', label='Equilibration Start')
    plt.axvline(x=plateau_results['mean_representative_time'], color='b', linestyle='--', label='Equilibration Point')
    plt.axhline(y=plateau_results['plateau_average'], 
            color='g', linestyle='--', 
            label=f'Plateau Average: {plateau_results["plateau_average"]:.2f} Å')

    # Add plateau region shading
    plt.fill_between(time[plateau_results['start_index']:],
                    plateau_results['plateau_average'] - plateau_results['plateau_std'],
                    plateau_results['plateau_average'] + plateau_results['plateau_std'],
                    color='g', alpha=0.2,
                    label=f'Std Dev: ±{plateau_results["plateau_std"]:.2f} Å')

    plt.xlabel("Time (ps)")
    plt.ylabel(r'RMSD ($\AA$)')
    plt.title(f'RMSD for {os.path.basename(pdb_name).replace(".pdb", "")}')
    plt.legend()
    # Make a new folder called 'rmsd_plots' to save the plots
    rmsd_plot_dir = os.path.join(output_dir, 'rmsd_plots')
    os.makedirs(rmsd_plot_dir, exist_ok=True)
    rmsd_plot_name = os.path.join(rmsd_plot_dir, f'{os.path.basename(pdb_name).replace(".pdb", "")}_mds.png')
    plt.savefig(rmsd_plot_name)
    plt.close()

    plateau_index = plateau_index + 1

    print(f"Equilibration analysis for {os.path.basename(pdb_name)}:")
    print(f"  Time: {plateau_time:.2f} ps = Model {plateau_index}")
    print(f"  Plateau average: {plateau_results['plateau_average']:.2f} Å")
    print(f"  Plateau std dev: {plateau_results['plateau_std']:.2f} Å")
    
    # Extract equilibration model PDB
    equilibration_pdb_name = os.path.join(output_dir, os.path.basename(pdb_name))
    
    # Using MDAnalysis to write the specific model
    u_eq = mda.Universe(mds_output_name or processed_filename)
    u_eq.trajectory[plateau_index]
    protein = u_eq.select_atoms("protein")
    
    with mda.Writer(equilibration_pdb_name, protein.n_atoms) as W:
        W.write(protein)
    
    # Close the second Universe
    u_eq.trajectory.close()
    del u_eq

    return processed_filename, mds_output_name, plateau_time, equilibration_pdb_name

def prep_receptors():
    """
    Process all PDB files in a specified directory
    
    Parameters:
    -----------
    directory : str
        Path to directory containing PDB files
    output_dir : str, optional
        Path to directory for saving outputs (default: 'prepared_structures')
    """
    # Set working directory
    pdb_dir = input("Enter the path to the directory containing the PDB files: ")
    # Check if the path exists
    while not os.path.exists(pdb_dir):
        pdb_dir = input(
            "That path does not appear to exist.\nPlease enter the path to "
            "the directory containing the PDB files: "
        )
    os.chdir(pdb_dir)

    # Create output directory
    output_dir = 'prepared_structures'
    os.makedirs(output_dir, exist_ok=True)
    
    # Ask user if a reference structure is present and if so, what is the name
    ref_structure = input("Is there a reference structure in the specified PDB directory? (y/n): ").lower().strip() == 'y'
    if ref_structure:
        ref_name = input("Enter the name of the reference structure: ")
        if not ref_name.endswith(".pdb"):
            ref_name += ".pdb"
        if ref_name not in os.listdir(pdb_dir):
            print(f"Reference structure {ref_name} not found in the directory.")
            ref_name = None
        # Add a copy of the reference structure to the output directory
        ref_path = os.path.join(pdb_dir, ref_name)
        shutil.copy(ref_path, os.path.join(output_dir, ref_name))
    else:
        ref_name = None

    # Get all PDB files in directory
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]

    # Remove the reference structure from the list
    if ref_structure:
        pdb_files.remove(ref_name)
    
    # Ask about saving trajectories
    print(f"\nFound {len(pdb_files)} PDB files in {pdb_dir}")
    perform_mds = input("Do you want to perform Molecular Dynamics Simulation on the structures? (y/n): ").lower().strip() == 'y'

    if not perform_mds:
        # Just run PDBFixer
        print("\nRunning PDBFixer on all PDB files...")
        results = []
        for pdb_file in pdb_files:
            # Skip the reference structure, if any
            if pdb_file.startswith("ref_"):
                print(f"Skipping reference structure: {pdb_file}")
                continue
            full_path = os.path.join(pdb_dir, pdb_file)
            try:
                result = run_pdbfixer(full_path, output_dir)
                results.append(result)
            except Exception as e:
                print(f"Failed to process {pdb_file} with PDBFixer. Error: {e}")
        print("\nAll PDB files have been processed with PDBFixer.")

        return results
    
    # If the user wants to perform Molecular Dynamics Simulation
    save_trajectories = input(f"Do you want to save the MDS trajectory files? (y/n): ").lower().strip() == 'y'
    mds_time = input("Enter the total simulation time in ns (10-100 ns): ")
    while not mds_time.isdigit() or int(mds_time) <= 0 or int(mds_time) > 100:
        mds_time = input("Please enter a valid time (10-100 ns): ")
    mds_time = int(mds_time)

    results = []
    for pdb_file in pdb_files:
        full_path = os.path.join(pdb_dir, pdb_file)
        try:
            processed_filename = run_pdbfixer(full_path, output_dir)
            result = run_mds(processed_filename, output_dir, mds_time)
            results.append(result)
            if not save_trajectories:
                os.remove(os.path.join(output_dir, pdb_file.replace(".pdb", "-mds.pdb")))
        except Exception as e:
            print(f"Failed to process {pdb_file} with MDS. Error: {e}")

    print("\nAll PDB files have been processed with Molecular Dynamics Simulation.")
    return results

if __name__ == "__main__":
    prep_receptors()