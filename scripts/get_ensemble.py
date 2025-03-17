import os
import json
import requests
import argparse
import csv
import re
import sys
from typing import Optional
import concurrent.futures
from typing import List, Optional, Tuple, Union, Dict, Any
from datetime import datetime
from collections import defaultdict

class PDBDownloader:
    """
    Class to search and download PDB structures from the RCSB PDB database
    based on user-defined criteria, with parallelization for improved performance.
    """

    def __init__(self, download_dir: str = "pdb_downloads", max_threads: int = 10):
        """
        Initialize the PDB Downloader with a download directory.

        Args:
            download_dir: Directory to save downloaded PDB files
            max_threads: Maximum number of concurrent download threads
        """
        self.download_dir = download_dir
        self.max_threads = max_threads
        self.search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        self.download_url = "https://files.rcsb.org/download"
        self.session = requests.Session()  # Using a session for connection pooling
        
        # Create download directory if it doesn't exist
        if not os.path.exists(self.download_dir):
            os.makedirs(self.download_dir)

    def build_query(self,
                    gene_name: Optional[Union[str, List[str]]] = None,
                    species: Optional[Union[str, List[str]]] = None,
                    terms: Optional[Union[str, List[str]]] = None,
                    mutations: int = 1,
                    ref_pdb_id: Optional[str] = None,
                    max_res: float = 2.5) -> Dict[str, Any]:
        """
        Build a query JSON for the RCSB Search API.

        Args:
            gene_name: Gene name to search for
            species: Scientific name(s) of the source organism(s). Can be a single string or a list of strings.
            mutations: Whether to include mutants (if False, only structures with mutation count < 1)
            ref_pdb_id: PDB ID of a reference structure for structural similarity search
            terms: Text search term(s) for full-text search. Can be a single string or a list of strings.
            max_res: Maximum refinement resolution (default 2.5Å)

        Returns:
            dict: Query JSON for the RCSB Search API
        """
        # Build the query nodes
        query_nodes = []
        
        # Add gene name filter if provided
        if gene_name:
            # Convert single gene name to list for consistent handling
            if isinstance(gene_name, str):
                gene_list = [gene_name]
            else:
                gene_list = gene_name

            # If we have multiple genes, use an "or" group
            if len(gene_list) > 1:
                gene_nodes = []
                for gene in gene_list:
                    gene_query = {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                            "operator": "exact_match",
                            "negation": False,
                            "value": gene
                        }
                    }
                    gene_nodes.append(gene_query)
                
                # Group all gene queries with "or" operator
                gene_group = {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": gene_nodes
                }
                query_nodes.append(gene_group)
            elif len(gene_list) == 1:
                # Only one gene, add it directly to the main query
                gene_query = {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "negation": False,
                        "value": gene_list[0]
                    }
                }
                query_nodes.append(gene_query)

        # Add full-text search terms if provided
        if terms:
            # Convert single term to list for consistent handling
            if isinstance(terms, str):
                terms_list = [terms]
            else:
                terms_list = terms
                
            # If we have multiple search terms, use an "or" group
            if len(terms_list) > 1:
                term_nodes = []
                for term in terms_list:
                    term_query = {
                        "type": "terminal",
                        "service": "full_text",
                        "parameters": {
                            "value": term
                        }
                    }
                    term_nodes.append(term_query)
                
                # Group all term queries with "or" operator
                term_group = {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": term_nodes
                }
                query_nodes.append(term_group)
            elif len(terms_list) == 1:
                # Only one term, add it directly to the main query
                term_query = {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": terms_list[0]
                    }
                }
                query_nodes.append(term_query)
                
        # Add species filter if provided
        if species:
            # Convert single species to list for consistent handling
            if isinstance(species, str):
                species_list = [species]
            else:
                species_list = species

            # If we have multiple species, use an "or" group
            if len(species_list) > 1:
                species_nodes = []
                for sp in species_list:
                    species_query = {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                            "operator": "exact_match",
                            "negation": False,
                            "value": sp
                        }
                    }
                    species_nodes.append(species_query)
                
                # Group all species queries with "or" operator
                species_group = {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": species_nodes
                }
                query_nodes.append(species_group)
            elif len(species_list) == 1:
                # Only one species, add it directly to the main query
                species_query = {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_scientific_name",
                        "operator": "exact_match",
                        "negation": False,
                        "value": species_list[0]
                    }
                }
                query_nodes.append(species_query)
            # If species is an empty list, we don't add any species filter

        # Add resolution filter
        res_query = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal",
                "negation": False,
                "value": max_res
            }
        }
        query_nodes.append(res_query)

        # Filter for non-mutants if mutations is False
        mutation_query = {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "entity_poly.rcsb_mutation_count",
                "operator": "less_or_equal",
                "negation": False,
                "value": mutations
            }
        }
        query_nodes.append(mutation_query)

        # Create the main group of criteria
        main_query_group = {
            "type": "group",
            "nodes": query_nodes,
            "logical_operator": "and"
        }
        
        # Start with base query
        query = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [main_query_group]
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 1000},
                "results_content_type": ["experimental"],
                "sort": [{"sort_by": "score", "direction": "desc"}],
                "scoring_strategy": "combined"
            }
        }

        # Parse reference PDB ID safely
        if ref_pdb_id:
            parts = ref_pdb_id.split(":")
            if len(parts) > 1:
                entry_id = parts[0]
                asym_id = parts[1]
            else:
                # Default to the full ID as both entry_id and asym_id if format is unexpected
                entry_id = ref_pdb_id
                asym_id = "A"  # Default chain ID
                
            similarity_query = {
                "type": "terminal",
                "service": "structure",
                "parameters": {
                    "operator": "strict_shape_match",
                    "target_search_space": "polymer_entity_instance",
                    "value": {
                        "entry_id": entry_id,
                        "asym_id": asym_id
                    }
                }
            }
            query["query"]["nodes"].append(similarity_query)
            
        return query

    def search_structures(self, query: dict) -> List[str]:
        """
        Execute search with the given query and return matching PDB IDs.

        Args:
            query: Query JSON for the RCSB Search API

        Returns:
            List[str]: List of PDB IDs matching the search criteria
        """
        try:
            headers = {
                "Content-Type": "application/json",
                "Accept": "application/json"
            }
            
            response = self.session.post(self.search_url, json=query, headers=headers)
            response.raise_for_status()  # Raise exception for HTTP errors
            
            result = response.json()
            
            # Extract PDB IDs from the response
            pdb_ids = [hit["identifier"] for hit in result.get("result_set", [])]
            
            return pdb_ids
        
        except requests.exceptions.RequestException as e:
            print(f"Error during search: {e}")
            print(f"Response content: {response.text if 'response' in locals() else 'No response'}")
            return []

    def save_pdb_ids_to_file(self, pdb_ids: List[str]) -> str:
        """
        Save the PDB IDs to a text file, one ID per line.

        Args:
            pdb_ids: List of PDB IDs to save

        Returns:
            str: Path to the saved file
        """
        # Create a timestamped filename in the download directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file_path = os.path.join(self.download_dir, f"pdb_ids_{timestamp}.txt")
        
        # Write the PDB IDs to the file, one per line
        with open(output_file_path, 'w') as f:
            for pdb_id in pdb_ids:
                f.write(f"{pdb_id}\n")
        
        return output_file_path

    def download_pdb(self, pdb_id: str, file_format: str = "pdb", max_retries: int = 3) -> Tuple[str, str]:
        """
        Download a PDB structure file with retry logic. If the initial download fails,
        try to download the CIF format instead.

        Args:
            pdb_id: PDB ID of the structure to download
            file_format: Format of the structure file (pdb, cif, etc.)
            max_retries: Maximum number of download attempts

        Returns:
            Tuple[str, str]: (PDB ID, Path to the downloaded file or empty string if failed)
        """
        # Determine file extension based on format
        primary_extension = file_format
        
        # Construct download URL and local file path
        local_path = os.path.join(self.download_dir, f"{pdb_id}.{primary_extension}")
        
        # Check if file already exists (to avoid redundant downloads)
        if os.path.exists(local_path) and os.path.getsize(local_path) > 0:
            print(f"File already exists: {local_path}")
            return pdb_id, local_path
        
        # Try primary format first
        primary_url = f"{self.download_url}/{pdb_id}.{primary_extension}"
        
        try:
            print(f"Downloading {pdb_id} from {primary_url}")
            response = self.session.get(primary_url, timeout=30)
            response.raise_for_status()
            
            # Save the file
            with open(local_path, "wb") as f:
                f.write(response.content)
            
            print(f"Successfully downloaded {pdb_id} to {local_path}")
            return pdb_id, local_path
        
        except requests.exceptions.RequestException as e:
            print(f"Error downloading {pdb_id} in {primary_extension} format: {e}")
            print(f"Trying CIF format instead...")
            
            # If primary format fails, try CIF format
            if primary_extension != "cif":
                cif_url = f"{self.download_url}/{pdb_id}.cif"
                cif_local_path = os.path.join(self.download_dir, f"{pdb_id}.cif")
                
                try:
                    print(f"Downloading {pdb_id} from {cif_url}")
                    response = self.session.get(cif_url, timeout=30)
                    response.raise_for_status()
                    
                    # Save the CIF file
                    with open(cif_local_path, "wb") as f:
                        f.write(response.content)
                    
                    print(f"Successfully downloaded {pdb_id} as CIF to {cif_local_path}")
                    return pdb_id, cif_local_path
                
                except requests.exceptions.RequestException as e:
                    print(f"Error downloading {pdb_id} in CIF format: {e}")
        
        # If both primary format and CIF fail, try alternative URLs as a last resort
        alternative_urls = [
            f"https://files.rcsb.org/view/{pdb_id}.{primary_extension}",  # Alternate URL
            f"https://models.rcsb.org/v1/{pdb_id}/model"  # Models API (for some newer entries)
        ]
        
        for attempt, url in enumerate(alternative_urls, start=1):
            try:
                print(f"Attempting alternate URL {attempt}/{len(alternative_urls)}: {url}")
                response = self.session.get(url, timeout=30)
                response.raise_for_status()
                
                # Save the file
                with open(local_path, "wb") as f:
                    f.write(response.content)
                
                print(f"Successfully downloaded {pdb_id} to {local_path}")
                return pdb_id, local_path
            
            except requests.exceptions.RequestException as e:
                print(f"Error downloading {pdb_id} from alternate URL: {e}")
                if attempt == len(alternative_urls):
                    print(f"Failed to download {pdb_id} after trying all alternatives")
                    return pdb_id, ""
                print(f"Trying next URL...")
        
        return pdb_id, ""

    def download_structures_parallel(self, pdb_ids: List[str], file_format: str = "pdb", max_retries: int = 3) -> Tuple[List[str], List[str]]:
        """
        Download multiple PDB structures in parallel.

        Args:
            pdb_ids: List of PDB IDs to download
            file_format: Format of the structure files (pdb, cif, etc.)
            max_retries: Maximum number of download attempts per structure

        Returns:
            Tuple[List[str], List[str]]: Tuple containing (successfully downloaded PDB IDs, failed PDB IDs)
        """
        downloaded_pdb_ids = []
        failed_pdb_ids = []
        total = len(pdb_ids)
        
        # Use ThreadPoolExecutor to download multiple structures in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_threads) as executor:
            # Start the downloads in parallel
            future_to_pdb = {
                executor.submit(self.download_pdb, pdb_id, file_format, max_retries): pdb_id 
                for pdb_id in pdb_ids
            }
            
            # Process the results as they complete
            completed = 0
            for future in concurrent.futures.as_completed(future_to_pdb):
                completed += 1
                pdb_id, file_path = future.result()
                
                print(f"Progress: [{completed}/{total}] - {pdb_id}")
                
                if file_path:
                    downloaded_pdb_ids.append(pdb_id)
                else:
                    failed_pdb_ids.append(pdb_id)
        
        # Summary of results
        if failed_pdb_ids:
            print(f"\nWARNING: Failed to download {len(failed_pdb_ids)} structures: {', '.join(failed_pdb_ids[:10])}" + 
                 (f" and {len(failed_pdb_ids) - 10} more" if len(failed_pdb_ids) > 10 else ""))
        
        print(f"\nSuccessfully downloaded {len(downloaded_pdb_ids)} out of {total} structures")
        
        return downloaded_pdb_ids, failed_pdb_ids

    def save_download_results_to_csv(self, downloaded_pdb_ids: List[str], failed_pdb_ids: List[str]) -> Optional[str]:
        """
        Save the download results to a CSV file with 'downloaded' and 'failed' columns if there are failed downloads.
        The file is automatically saved to the download directory.

        Args:
            downloaded_pdb_ids: List of successfully downloaded PDB IDs
            failed_pdb_ids: List of PDB IDs that failed to download

        Returns:
            Optional[str]: Path to the created CSV file, or None if no CSV was created
        """
        if not failed_pdb_ids:
            print("No failed downloads, CSV will not be saved.")
            return None

        # Create a timestamped filename in the download directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file_path = os.path.join(self.download_dir, f"download_results_{timestamp}.csv")

        # Determine the maximum length needed for rows
        max_rows = max(len(downloaded_pdb_ids), len(failed_pdb_ids))
        
        with open(output_file_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write header
            writer.writerow(['downloaded', 'failed'])
            
            # Write data rows
            for i in range(max_rows):
                row = []
                # Add downloaded PDB ID if available for this row
                if i < len(downloaded_pdb_ids):
                    row.append(downloaded_pdb_ids[i])
                else:
                    row.append('')
                
                # Add failed PDB ID if available for this row
                if i < len(failed_pdb_ids):
                    row.append(failed_pdb_ids[i])
                else:
                    row.append('')
                
                writer.writerow(row)
        
        print(f"\nSaved download results to: {output_file_path}")
        return output_file_path

from collections import defaultdict
import os
from Bio.PDB import MMCIFParser, MMCIF2Dict, PDBIO

def count_mutations(structure_file):
    """
    Count engineered mutations in a PDB or CIF file.
    
    For PDB files: Analyzes SEQADV records for engineered mutations.
    For CIF files: Analyzes _struct_ref_seq_dif records for engineered mutations.
    
    Parameters:
    -----------
    structure_file : str
        Path to the PDB or CIF file
        
    Returns:
    --------
    dict
        Dictionary with chain IDs as keys and number of mutations as values
    """
    mutations = defaultdict(int)
    
    try:
        file_extension = os.path.splitext(structure_file)[1].lower()
        
        if file_extension == '.pdb':
            # Process PDB file
            with open(structure_file, 'r') as f:
                lines = f.readlines()
                
            # Count engineered mutations from SEQADV records
            for line in lines:
                if line.startswith("SEQADV") and "ENGINEERED MUTATION" in line:
                    chain_id = line[16:17].strip()  # Chain ID is at position 16
                    mutations[chain_id] += 1
                    
        elif file_extension == '.cif':
            # Process CIF file using BioPython
            cif_dict = MMCIF2Dict.MMCIF2Dict(structure_file)
            
            # Check if the necessary category exists
            if '_struct_ref_seq_dif.details' in cif_dict:
                details = cif_dict['_struct_ref_seq_dif.details']
                
                # Get the chains from the correct field
                # The chain ID in mmCIF is typically stored in _struct_ref_seq_dif.pdbx_auth_seq_id
                # or _struct_ref_seq_dif.pdbx_pdb_strand_id depending on the file
                if '_struct_ref_seq_dif.pdbx_pdb_strand_id' in cif_dict:
                    chains = cif_dict['_struct_ref_seq_dif.pdbx_pdb_strand_id']
                elif '_struct_ref_seq_dif.chain_id' in cif_dict:
                    chains = cif_dict['_struct_ref_seq_dif.chain_id']
                else:
                    # Fallback option if standard fields aren't available
                    print(f"Warning: Chain ID field not found in {structure_file}")
                    # Print available fields for debugging
                    ref_seq_fields = [k for k in cif_dict.keys() if k.startswith('_struct_ref_seq_dif')]
                    if ref_seq_fields:
                        print(f"Available fields: {ref_seq_fields}")
                    return dict(mutations)
                
                # Count engineered mutations
                for i, detail in enumerate(details):
                    if 'ENGINEERED MUTATION' in detail.upper():
                        chain_id = chains[i].strip()
                        mutations[chain_id] += 1
        else:
            print(f"Unsupported file format: {file_extension}")
            
    except Exception as e:
        print(f"Error parsing {structure_file}: {e}")
    
    return dict(mutations)

def id_pdb_mutants(structure_directory, mutation_threshold=1):
    """
    Screen a directory of PDB/CIF files and identify those with mutation counts exceeding a threshold.
    
    Parameters:
    -----------
    structure_directory : str
        Path to directory containing PDB and/or CIF files
    mutation_threshold : int
        Threshold for mutation count, files with more mutations than this will be reported
        
    Returns:
    --------
    dict
        Dictionary with structure file names (without extension) as keys and 
        dictionaries of mutation counts per chain as values
    """
    results = {}
    
    # Get all PDB and CIF files in the directory
    structure_files = [f for f in os.listdir(structure_directory) 
                      if f.endswith('.pdb') or f.endswith('.cif')]
    
    if not structure_files:
        print(f"No PDB or CIF files found in {structure_directory}")
        return results
    
    print(f"Analyzing {len(structure_files)} structure files for mutations...")
    
    for structure_file in structure_files:
        structure_id = os.path.splitext(structure_file)[0]
        file_path = os.path.join(structure_directory, structure_file)
        
        # Count mutations
        mutation_counts = count_mutations(file_path)
        
        # Check if any chain exceeds the threshold
        if any(count > mutation_threshold for count in mutation_counts.values()):
            results[structure_id] = mutation_counts
    
    return results

def convert_cif_to_pdb(cif_file: str, output_dir: str) -> Optional[str]:
    """
    Convert a structure in mmCIF format to PDB format using custom parsing.
    Function based on pdb_fromcif.py from the pdb-tools package.
    
    Args:
        cif_file: Path to the CIF file
        output_dir: Directory to save the converted PDB file
    
    Returns:
        Path to the converted PDB file or None if conversion fails
    """
    try:
        # Create output filename
        pdb_file_path = os.path.join(output_dir, os.path.splitext(os.path.basename(cif_file))[0] + ".pdb")
        
        with open(cif_file, 'r') as f_in, open(pdb_file_path, 'w') as f_out:
            # Format string for PDB ATOM/HETATM records
            _a = "{:6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}"
            _a += "{:6.2f}{:6.2f}      {:<4s}{:<2s}{:2s}\n"

            in_section, read_atom = False, False
            label_pos = 0
            labels = {}
            empty = set(('.', '?'))
            prev_model = None
            atom_num = 0
            serial = 0
            model_data = []  # store atom data to account for multi-model files
            
            for line in f_in:
                if line.startswith('loop_'):  # start of section
                    in_section = True

                elif line.startswith('#'):  # end of section
                    in_section = False
                    read_atom = False

                elif in_section and line.startswith('_atom_site.'):  # ATOM/HETATM
                    read_atom = True
                    labels[line.strip()] = label_pos
                    label_pos += 1

                elif read_atom and line.startswith(('ATOM', 'HETATM')):  # convert
                    fields = re.findall(r'[^"\s]\S*|".+?"', line)  # find enclosed ''

                    # Pick fields, giving preference to auth to match PDBs
                    model_no = fields[labels.get('_atom_site.pdbx_PDB_model_num')]
                    if prev_model != model_no:  # first line will trigger
                        prev_model = model_no
                        model_data.append([])
                        serial = 0

                    record = fields[labels.get('_atom_site.group_PDB')]
                    serial += 1

                    fid = labels.get('_atom_site.auth_atom_id')
                    if fid is None:
                        fid = labels.get('_atom_site.label_atom_id')
                    atname = fields[fid]

                    element = fields[labels.get('_atom_site.type_symbol')]
                    if element in empty:
                        element = ' '

                    # handle atom name
                    if atname[0] == '"' and atname[-1] == '"':
                        atname = atname[1:-1]

                    if len(atname) < 4 and atname[0].isalpha() and len(element) < 2:
                        atname = ' ' + atname  # pad

                    altloc = fields[labels.get('_atom_site.label_alt_id')]
                    if altloc in empty:
                        altloc = ' '

                    fid = labels.get('_atom_site.auth_comp_id')
                    if fid is None:
                        fid = labels.get('_atom_site.label_comp_id')
                    resname = fields[fid]

                    fid = labels.get('_atom_site.auth_asym_id')
                    if fid is None:
                        fid = labels.get('_atom_site.label_asym_id')
                    chainid = fields[fid]

                    fid = labels.get('_atom_site.auth_seq_id')
                    if fid is None:
                        fid = labels.get('_atom_site.label_seq_id')
                    resnum = int(fields[fid])

                    icode = fields[labels.get('_atom_site.pdbx_PDB_ins_code')]
                    if icode in empty:
                        icode = ' '

                    x = float(fields[labels.get('_atom_site.Cartn_x')])
                    y = float(fields[labels.get('_atom_site.Cartn_y')])
                    z = float(fields[labels.get('_atom_site.Cartn_z')])
                    occ = float(fields[labels.get('_atom_site.occupancy')])
                    bfactor = float(fields[labels.get('_atom_site.B_iso_or_equiv')])

                    charge = fields[labels.get('_atom_site.pdbx_formal_charge')]
                    try:
                        charge = charge
                    except ValueError:
                        charge = '  '

                    segid = chainid

                    atom_line = _a.format(record, serial, atname, altloc, resname,
                                          chainid, resnum, icode, x, y, z, occ, bfactor,
                                          segid, element, charge)

                    atom_num += 1

                    # Check if structure is too large
                    if atom_num > 99999:
                        raise ValueError(f"Number of atoms exceeds PDB format limit: {atom_num}")
                    elif len(chainid) > 1:
                        raise ValueError(f"Chain ID is too large: {chainid}")
                    elif resnum > 9999:
                        raise ValueError(f"Too many residues ({resnum}) in chain {chainid}")

                    model_data[-1].append(atom_line)

            # Write PDB data to output file
            is_ensemble = len(model_data) > 1
            if is_ensemble:
                for model_no, model in enumerate(model_data, start=1):
                    f_out.write("MODEL {:>5d}\n".format(model_no))
                    for line in model:
                        f_out.write(line)
                    f_out.write('ENDMDL\n')
            else:
                for line in model_data[0]:
                    f_out.write(line)

            f_out.write("{:<80s}\n".format("END"))
            
        return pdb_file_path
        
    except Exception as e:
        print(f"Error converting {cif_file} to PDB: {e}")
        return None

def main():
    """Command-line interface for the PDB Downloader."""
    parser = argparse.ArgumentParser(description="Download PDB structures based on search criteria")
    
    parser.add_argument("--gene", type=str, nargs="+", help="Gene name(s). Multiple genes can be provided.")
    parser.add_argument("--species", type=str, nargs="+", help="Scientific name(s) of the source organism(s). Multiple species can be provided.")
    parser.add_argument("--terms", type=str, nargs="+", help="Full-text search terms. Multiple terms can be provided.")
    parser.add_argument("--mutations", type=int, default=1, help="Maximum number of allowed mutations (default: 1)")
    parser.add_argument("--ref", type=str, help="Reference PDB ID for structural similarity (include chain ID, e.g., 1A2B:A)")
    parser.add_argument("--res", type=float, default=2.5, help="Maximum resolution (default: 2.5Å)")
    parser.add_argument("--format", type=str, default="pdb", choices=["pdb", "cif", "xml"], 
                       help="File format for download (default: pdb)")
    parser.add_argument("--output-dir", type=str, default="pdb_downloads", 
                       help="Directory to save downloaded files (default: pdb_downloads)")
    parser.add_argument("--limit", type=int, default=1000, 
                       help="Maximum number of structures to download (default: 1000)")
    parser.add_argument("--threads", type=int, default=10,
                       help="Number of parallel download threads (default: 10)")
    parser.add_argument("--skip-existing", action="store_true",
                       help="Skip downloading files that already exist")
    parser.add_argument("--yes", "-y", action="store_true",
                       help="Skip confirmation prompt and proceed with downloads")
    parser.add_argument("--ids-only", action="store_true",
                       help="Output PDB IDs to a file without downloading structures")
    # Add output file argument for mutation analysis results
    parser.add_argument("--output-mutations", type=str, 
                       help="Output file for mutation analysis results (default: mutations_TIMESTAMP.csv)")
    
    args = parser.parse_args()
    
    # Initialize the downloader with thread count
    downloader = PDBDownloader(download_dir=args.output_dir, max_threads=args.threads)
    
    # Build the query
    query = downloader.build_query(
        gene_name=args.gene,
        species=args.species,
        terms=args.terms,
        mutations=args.mutations,
        ref_pdb_id=args.ref,
        max_res=args.res
    )
    
    # Update pagination limit
    query["request_options"]["paginate"]["rows"] = args.limit
    
    # Print the query for debugging
    print("Search query:")
    print(json.dumps(query, indent=2))
    
    # Always save the query to a file in the output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create filename components
    filename_parts = ["query"]
    
    # Add gene information if provided
    if args.gene:
        gene_str = "_".join([g.replace(' ', '-') for g in args.gene])
        filename_parts.append(f"gene-{gene_str}")
    
    # Add species information if provided
    if args.species:
        species_str = "_".join([sp.replace(' ', '-') for sp in args.species])
        filename_parts.append(f"species-{species_str}")
    
    # Add terms information if provided
    if args.terms:
        terms_str = "_".join([t.replace(' ', '-') for t in args.terms])
        filename_parts.append(f"terms-{terms_str}")
    
    # Add reference if provided
    if args.ref:
        ref_id = args.ref.split(":")[0]
        filename_parts.append(f"ref-{ref_id}")
    
    # Add timestamp
    filename_parts.append(timestamp)
    
    # Make sure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
        
    # Join all parts with underscores
    query_file = "_".join(filename_parts) + ".json"

    try:
        query_file_path = os.path.join(args.output_dir, query_file)
        with open(query_file_path, 'w') as f:
            json.dump(query, f, indent=2)
        print(f"Saved search query to {query_file_path}")
    except IOError as e:
        print(f"Error saving query to file: {e}")

    # Execute the search
    print("\nSearching for structures...")
    pdb_ids = downloader.search_structures(query)
    
    if not pdb_ids:
        print("No structures found matching the criteria.")
        return
    
    print(f"\nFound {len(pdb_ids)} matching structures: {', '.join(pdb_ids[:10])}" + 
        (f" and {len(pdb_ids) - 10} more" if len(pdb_ids) > 10 else ""))
    
    # If --ids-only flag is set, save the IDs to a file and exit
    if args.ids_only:
        output_file_path = downloader.save_pdb_ids_to_file(pdb_ids)
        print(f"PDB IDs saved to {output_file_path}")
        return
    
    # Ask user to confirm download unless --yes flag is set
    if not args.yes:
        response = input(f"\nDo you want to download {len(pdb_ids)} structures? (y/n): ")
        if response.lower() != 'y':
            print("Download canceled.")
            return
    
    # Download the structures with parallel processing
    print(f"\nDownloading structures...")
    downloaded_pdb_ids, failed_pdb_ids = downloader.download_structures_parallel(
        pdb_ids, file_format=args.format, max_retries=3
    )
    
    # Analyze mutations after download using the mutations threshold from args, unless XML format
    if args.format == "xml":
        print("Skipping mutation analysis for XML format as it is not supported.")
        return
    print("\nAnalyzing structures for mutations...")
    
    # Create a consolidated results file
    results_file = args.output_mutations if args.output_mutations else f"pdb_results_{timestamp}.csv"
    results_file_path = os.path.join(args.output_dir, results_file)
    
    # Create a dictionary to store all PDB information
    pdb_info = {}
    
    # Initialize with all PDB IDs from the search
    for pdb_id in pdb_ids:
        pdb_info[pdb_id] = {
            "downloaded": pdb_id in downloaded_pdb_ids,
            "mutations": {},
            "highest_mutation_count": 0,
            "over_threshold": False,
            "action": "Kept"
        }
    
    # Run the mutation analysis using the same threshold as the download query
    mutation_results = id_pdb_mutants(args.output_dir, args.mutations)
    
    # Process mutation results
    structures_over_threshold = []
    
    for structure_id, chain_mutations in mutation_results.items():
        if structure_id in pdb_info:
            pdb_info[structure_id]["mutations"] = chain_mutations
            highest_count = max(chain_mutations.values()) if chain_mutations else 0
            pdb_info[structure_id]["highest_mutation_count"] = highest_count
            pdb_info[structure_id]["over_threshold"] = highest_count > args.mutations
            
            if pdb_info[structure_id]["over_threshold"]:
                structures_over_threshold.append(structure_id)
    
    # If structures with mutations exceeding threshold are found, ask user if they want to keep them
    keep_structures = True
    if structures_over_threshold:
        print(f"\nFound {len(structures_over_threshold)} structures with mutations exceeding threshold {args.mutations}:")
        
        # Display structures with excessive mutations
        for pdb_id in structures_over_threshold:
            chain_mutations = pdb_info[pdb_id]["mutations"]
            for chain, count in chain_mutations.items():
                if count > args.mutations:
                    print(f"  {pdb_id} Chain {chain}: {count} mutations")
        
        # Ask user if they want to keep these structures
        if not args.yes:  # If not in auto-yes mode
            response = input(f"\nDo you want to keep these {len(structures_over_threshold)} structures with excessive mutations? (y/n): ")
            keep_structures = response.lower() == 'y'
        
        # Process user's choice
        if not keep_structures:
            for pdb_id in structures_over_threshold:
                pdb_info[pdb_id]["action"] = "Removed"
                
                # Remove the file if the user doesn't want to keep it
                try:
                    file_path = os.path.join(args.output_dir, f"{pdb_id}.{args.format}")
                    if os.path.exists(file_path):
                        os.remove(file_path)
                        print(f"Removed {file_path}")
                except Exception as e:
                    print(f"Error removing {pdb_id}: {e}")
                    pdb_info[pdb_id]["action"] = "Failed to remove"
            
            print(f"\nRemoved {len(structures_over_threshold)} structures with excessive mutations")
    else:
        print(f"No structures found with more than {args.mutations} mutations.")
    
    # Write consolidated results to CSV
    with open(results_file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header
        writer.writerow(['PDB_ID', 'Downloaded', 'Mutations_by_Chain', 
                        'Mutation_Threshold', 'Over_Threshold', 'Action'])
        
        # Write data rows
        for pdb_id, info in pdb_info.items():
            mutations_str = "; ".join([f"{chain}:{count}" for chain, count in info["mutations"].items()]) if info["mutations"] else "NA"
            writer.writerow([
                pdb_id,
                "Yes" if info["downloaded"] else "No",
                mutations_str,
                args.mutations,
                "Yes" if info["over_threshold"] else "No",
                info["action"]
            ])
    
    print(f"\nDownload and mutation analysis results saved to {results_file_path}")

    # If there are any CIF files in the output directory, ask user if they want to convert them to PDB format
    if any(file.endswith('.cif') for file in os.listdir(args.output_dir)):
        response = input(f"\nDo you want to convert CIF files to PDB format? (y/n): ")
        if response.lower() == 'y':
            print("Converting CIF files to PDB format...")
            for cif_file in os.listdir(args.output_dir):
                if cif_file.endswith('.cif'):
                    cif_path = os.path.join(args.output_dir, cif_file)
                    pdb_path = convert_cif_to_pdb(cif_path, args.output_dir)
                    if pdb_path:
                        print(f"Converted {cif_file} to {pdb_path}")
                    else:
                        print(f"Failed to convert {cif_file}")

if __name__ == "__main__":
    main()