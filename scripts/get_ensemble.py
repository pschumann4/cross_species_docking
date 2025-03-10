import os
import json
import requests
import argparse
import csv
import time
import concurrent.futures
from typing import List, Optional, Tuple, Union, Dict, Any
from datetime import datetime

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
                    include_mutants: bool = False,
                    ref_pdb_id: Optional[str] = None,
                    max_res: float = 2.5) -> Dict[str, Any]:
        """
        Build a query JSON for the RCSB Search API.

        Args:
            gene_name: Gene name to search for
            species: Scientific name(s) of the source organism(s). Can be a single string or a list of strings.
            include_mutants: Whether to include mutants (if False, only structures with mutation count = 0)
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

        # Filter for non-mutants if include_mutants is False
        if not include_mutants:
            mutation_query = {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "entity_poly.rcsb_mutation_count",
                    "operator": "equals",
                    "negation": False,
                    "value": 0
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

        # Add structural similarity search if ref_pdb_id is provided
        if ref_pdb_id:
            similarity_query = {
                "type": "terminal",
                "service": "structure",
                "parameters": {
                    "operator": "strict_shape_match",
                    "target_search_space": "assembly",
                    "value": {
                        "entry_id": ref_pdb_id,
                        "assembly_id": "1"
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

    def download_pdb(self, pdb_id: str, file_format: str = "pdb", max_retries: int = 3) -> Tuple[str, str]:
        """
        Download a PDB structure file with retry logic.

        Args:
            pdb_id: PDB ID of the structure to download
            file_format: Format of the structure file (pdb, cif, etc.)
            max_retries: Maximum number of download attempts

        Returns:
            Tuple[str, str]: (PDB ID, Path to the downloaded file or empty string if failed)
        """
        # Determine file extension based on format
        extension = file_format
        
        # Construct download URL and local file path
        local_path = os.path.join(self.download_dir, f"{pdb_id}.{extension}")
        
        # Check if file already exists (to avoid redundant downloads)
        if os.path.exists(local_path) and os.path.getsize(local_path) > 0:
            print(f"File already exists: {local_path}")
            return pdb_id, local_path
        
        # Try alternative URLs if initial one fails
        urls_to_try = [
            f"{self.download_url}/{pdb_id}.{extension}",  # Standard URL
            f"https://files.rcsb.org/view/{pdb_id}.{extension}",  # Alternate URL
            f"https://models.rcsb.org/v1/{pdb_id}/model"  # Models API (for some newer entries)
        ]
        
        for attempt in range(max_retries):
            try:
                # Use a different URL for each retry
                current_url = urls_to_try[min(attempt, len(urls_to_try) - 1)]
                
                print(f"Downloading {pdb_id} from {current_url}")
                response = self.session.get(current_url, timeout=30)
                response.raise_for_status()
                
                # Save the file
                with open(local_path, "wb") as f:
                    f.write(response.content)
                
                print(f"Successfully downloaded {pdb_id} to {local_path}")
                return pdb_id, local_path
            
            except requests.exceptions.RequestException as e:
                print(f"Error downloading {pdb_id} (attempt {attempt+1}/{max_retries}): {e}")
                if attempt == max_retries - 1:
                    print(f"Failed to download {pdb_id} after {max_retries} attempts")
                    return pdb_id, ""
                print(f"Retrying in 2 seconds...")
                time.sleep(2)  # Wait before retry
        
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


def main():
    """Command-line interface for the PDB Downloader."""
    parser = argparse.ArgumentParser(description="Download PDB structures based on search criteria")
    
    parser.add_argument("--gene", type=str, nargs="+", help="Gene name(s). Multiple genes can be provided.")
    parser.add_argument("--species", type=str, nargs="+", help="Scientific name(s) of the source organism(s). Multiple species can be provided.")
    parser.add_argument("--terms", type=str, nargs="+", help="Full-text search terms. Multiple terms can be provided (combined with OR logic).")
    parser.add_argument("--include-mutants", action="store_true", help="Include structures with mutations")
    parser.add_argument("--ref", type=str, help="Reference PDB ID for structural similarity")
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
    
    args = parser.parse_args()
    
    # Initialize the downloader with thread count
    downloader = PDBDownloader(download_dir=args.output_dir, max_threads=args.threads)
    
    # Build the query
    query = downloader.build_query(
        gene_name=args.gene,
        species=args.species,
        terms=args.terms,
        include_mutants=args.include_mutants,
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
        filename_parts.append(f"ref-{args.ref}")
    
    # Add timestamp
    filename_parts.append(timestamp)
    
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
    
    print(f"\nDownloaded {len(downloaded_pdb_ids)} structures to {args.output_dir}")
    
    # Save results to CSV if there were any failed downloads
    downloader.save_download_results_to_csv(downloaded_pdb_ids, failed_pdb_ids)


if __name__ == "__main__":
    main()