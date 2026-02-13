def run_tmalign():    
    """
    Calculate TM-scores across model directory using TM-align executable.

    Methodology:
    1. Load and preprocess data - identify PDB files in directory
    2. Identify reference model (prefixed with "ref_")
    3. Run TM-align for each test model against the reference
    4. Extract TM-scores (both normalizations) and alignment statistics
    5. Output results to CSV
    """
    import os
    import subprocess
    import pandas as pd
    import re

    # Prompt user for the PDB directory
    pdb_dir = input("Enter the directory containing the PDB models: ").strip()

    # Check that the directory exists
    while not os.path.exists(pdb_dir):
        print("The directory does not exist.")
        pdb_dir = input("Enter the directory containing the PDB models: ").strip()

    # Check for the presence of a reference PDB file (starting with "ref_")
    ref_pdb = None
    pdb_files = [i for i in os.listdir(pdb_dir) if i.endswith(".pdb")]
    
    for file in pdb_files:
        if file.startswith("ref_"):
            ref_pdb = file
            print(f"Found reference structure: {ref_pdb}")
            break
    
    while ref_pdb is None:
        # Ask for the filename of the reference PDB file
        ref_pdb = input("Please enter the filename of the reference PDB file (in the specified directory): ").strip()

        # Check if the file exists in the directory
        if not os.path.exists(os.path.join(pdb_dir, ref_pdb)):
            print(f"The reference PDB file '{ref_pdb}' does not exist in {pdb_dir}.")
            ref_pdb = None
        elif not ref_pdb.endswith(".pdb"):
            print("The reference file must be a PDB file (.pdb extension).")
            ref_pdb = None

    # Run TM-align for each test model against the reference
    results = []
    ref_pdb_path = os.path.join(pdb_dir, ref_pdb)
    
    print(f"\nRunning TM-align against reference: {ref_pdb}")
    print("-" * 60)
    
    for file in pdb_files:
        if file.endswith(".pdb") and file != ref_pdb:
            test_pdb = os.path.join(pdb_dir, file)
            
            # Run TM-align with the test structure as query (first arg) and reference as target (second arg)
            command = f'TMalign "{test_pdb}" "{ref_pdb_path}"'
            
            try:
                process = subprocess.Popen(
                    command, 
                    shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.PIPE,
                    cwd=pdb_dir  # Run in the PDB directory
                )
                stdout, stderr = process.communicate()
                
                if process.returncode != 0:
                    print(f"Warning: TM-align returned error for {file}")
                    print(f"Error: {stderr.decode()}")
                    continue
                
                # Parse TM-align output
                output = stdout.decode()
                
                # Extract key metrics from TM-align output
                tm_score_chain1 = None  # TM-score normalized by length of chain1 (test)
                tm_score_chain2 = None  # TM-score normalized by length of chain2 (reference)
                rmsd = None
                aligned_length = None
                seq_id = None
                chain1_length = None
                chain2_length = None
                
                for line in output.splitlines():
                    # Length of structures
                    if line.startswith("Length of Chain_1:"):
                        match = re.search(r'Length of Chain_1:\s*(\d+)', line)
                        if match:
                            chain1_length = int(match.group(1))
                    
                    elif line.startswith("Length of Chain_2:"):
                        match = re.search(r'Length of Chain_2:\s*(\d+)', line)
                        if match:
                            chain2_length = int(match.group(1))
                    
                    # Aligned length, RMSD, and Seq_ID are on the same line
                    elif line.startswith("Aligned length="):
                        # Extract all three values from this line
                        # Format: "Aligned length= 242, RMSD=   1.75, Seq_ID=n_identical/n_aligned= 0.921"
                        aligned_match = re.search(r'Aligned length=\s*(\d+)', line)
                        rmsd_match = re.search(r'RMSD=\s*([\d.]+)', line)
                        seqid_match = re.search(r'Seq_ID=n_identical/n_aligned=\s*([\d.]+)', line)
                        
                        if aligned_match:
                            aligned_length = int(aligned_match.group(1))
                        if rmsd_match:
                            rmsd = float(rmsd_match.group(1))
                        if seqid_match:
                            seq_id = float(seqid_match.group(1))
                    
                    # TM-score normalized by length of Chain_1 (test model)
                    elif line.startswith("TM-score=") and "Chain_1" in line:
                        # Format: "TM-score= 0.87116 (if normalized by length of Chain_1, i.e., LN=259, d0=5.95)"
                        match = re.search(r'TM-score=\s*([\d.]+)', line)
                        if match:
                            tm_score_chain1 = float(match.group(1))
                    
                    # TM-score normalized by length of Chain_2 (reference)
                    elif line.startswith("TM-score=") and "Chain_2" in line:
                        # Format: "TM-score= 0.92513 (if normalized by length of Chain_2, i.e., LN=243, d0=5.78)"
                        match = re.search(r'TM-score=\s*([\d.]+)', line)
                        if match:
                            tm_score_chain2 = float(match.group(1))
                
                results.append({
                    "Test_Model": file,
                    "Chain1_Length": chain1_length,
                    "Chain2_Length": chain2_length,
                    "TM-score_norm_test": tm_score_chain1,  # Normalized by test length
                    "TM-score_norm_ref": tm_score_chain2,   # Normalized by reference length
                    "RMSD": rmsd,
                    "Aligned_Length": aligned_length,
                    "Seq_Identity": seq_id
                })
                
                print(f"Processed: {file:40s} TM-score(ref): {tm_score_chain2:.5f}" if tm_score_chain2 else f"Processed: {file:40s} [parsing failed]")
                
            except Exception as e:
                print(f"Error processing {file}: {str(e)}")
                continue

    # Save results to a CSV file
    if results:
        results_df = pd.DataFrame(results)
        output_file = os.path.join(pdb_dir, "tm_scores.csv")
        results_df.to_csv(output_file, index=False)
        
        print("\n" + "=" * 60)
        print(f"TM-align analysis complete!")
        print(f"Results saved to: {output_file}")
        print(f"Total models analyzed: {len(results)}")
        print("=" * 60)
        
        # Print summary statistics
        if "TM-score_norm_ref" in results_df.columns:
            valid_scores = results_df["TM-score_norm_ref"].dropna()
            if len(valid_scores) > 0:
                print(f"\nTM-score summary (normalized by reference):")
                print(f"  Mean:   {valid_scores.mean():.5f}")
                print(f"  Median: {valid_scores.median():.5f}")
                print(f"  Min:    {valid_scores.min():.5f}")
                print(f"  Max:    {valid_scores.max():.5f}")
                print(f"  Std:    {valid_scores.std():.5f}")
    else:
        print("No results were generated. Please check your input files and TM-align installation.")

if __name__ == "__main__":
    run_tmalign()