import requests
import csv
import time
from time import sleep
from typing import List, Dict, Tuple, Optional

def search_pdb_by_smiles(smiles: str) -> List[str]:
    """
    Search RCSB PDB for a given SMILES string and return PDB IDs with exact matches.
    
    Args:
        smiles (str): SMILES string of the chemical compound
        
    Returns:
        list: List of PDB IDs with exact matches (score = 1)
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    query = {
        "query": {
            "type": "terminal",
            "service": "chemical",
            "parameters": {
                "value": smiles,
                "type": "descriptor",
                "descriptor_type": "SMILES",
                "match_type": "graph-strict"
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 1000
            },
            "results_content_type": [
                "experimental"
            ],
            "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
            ],
            "scoring_strategy": "combined"
        }
    }
    
    try:
        response = requests.post(url, json=query)
        response.raise_for_status()
        results = response.json()
        
        exact_matches = [
            hit['identifier'] 
            for hit in results.get('result_set', [])
            if hit.get('score', 0) == 1
        ]
        
        return exact_matches
        
    except requests.exceptions.RequestException as e:
        print(f"Error searching PDB: {e}")
        return []

def get_pdb_metadata(pdb_id: str) -> Optional[Dict]:
    """
    Retrieve metadata for a PDB entry using the RCSB Data API.
    
    This fetches:
    - Experimental method
    - Resolution (for X-ray structures)
    - Primary citation DOI
    - Structure title
    - Polymer entity information (for organism source)
    
    Args:
        pdb_id (str): PDB identifier
        
    Returns:
        dict: Metadata dictionary or None if request fails
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        # Extract basic metadata
        metadata = {
            'pdb_id': pdb_id,
            'title': data.get('struct', {}).get('title', 'N/A'),
            'method': 'N/A',
            'resolution': 'N/A',
            'db_link': f'https://www.rcsb.org/structure/{pdb_id}'
        }
        
        # Get experimental method
        exptl = data.get('exptl', [])
        if exptl and len(exptl) > 0:
            metadata['method'] = exptl[0].get('method', 'N/A')
        
        # Get resolution (primarily for X-ray structures)
        refine = data.get('refine', [])
        if refine and len(refine) > 0:
            resolution = refine[0].get('ls_d_res_high')
            if resolution:
                metadata['resolution'] = f"{resolution:.2f}"
        
        # Alternative resolution from rcsb_entry_info
        if metadata['resolution'] == 'N/A':
            resolution_combined = data.get('rcsb_entry_info', {}).get('resolution_combined')
            if resolution_combined and len(resolution_combined) > 0:
                metadata['resolution'] = f"{resolution_combined[0]:.2f}"
        
        return metadata
        
    except requests.exceptions.RequestException as e:
        print(f"Error fetching metadata for {pdb_id}: {e}")
        return None

def get_polymer_entities(pdb_id: str) -> List[Dict]:
    """
    Retrieve polymer entity information using GraphQL API.
    
    Uses the standard RCSB structure query format to extract:
    - Molecule name
    - Source organism
    - Sequence length
    - Mutation information
    
    Args:
        pdb_id (str): PDB identifier
        
    Returns:
        list: List of polymer entity dictionaries
    """
    url = "https://data.rcsb.org/graphql"
    
    # Using the standard RCSB structure query format
    query = """
    query structure($id: String!) {
      entry(entry_id: $id) {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            entry_id
            entity_id
            auth_asym_ids
          }
          rcsb_polymer_entity {
            pdbx_description
          }
          entity_poly {
            type
            rcsb_entity_polymer_type
            pdbx_seq_one_letter_code_can
            rcsb_sample_sequence_length
            rcsb_mutation_count
          }
          rcsb_entity_source_organism {
            scientific_name
            ncbi_scientific_name
          }
          rcsb_entity_host_organism {
            ncbi_scientific_name
          }
        }
      }
    }
    """
    
    try:
        response = requests.post(
            url,
            json={'query': query, 'variables': {'id': pdb_id}},
            timeout=30
        )
        response.raise_for_status()
        
        # Parse JSON response
        try:
            data = response.json()
        except ValueError as e:
            print(f"    Failed to parse JSON response for {pdb_id}: {e}")
            return []
        
        # Check for GraphQL errors
        if 'errors' in data:
            print(f"    GraphQL errors for {pdb_id}: {data['errors']}")
            return []
        
        # Check if response has expected structure
        if not data or 'data' not in data:
            print(f"    Empty or malformed response from GraphQL API for {pdb_id}")
            return []
        
        if not data['data'] or 'entry' not in data['data']:
            print(f"    No entry data found for {pdb_id}")
            return []
        
        if not data['data']['entry']:
            print(f"    Entry is null for {pdb_id} (may not exist or be obsolete)")
            return []
        
        entities = []
        polymer_entities = data['data']['entry'].get('polymer_entities', [])
        
        if not polymer_entities:
            print(f"    No polymer entities found for {pdb_id}")
            return []
        
        for entity in polymer_entities:
            # Get chain identifiers
            chains = 'N/A'
            if entity.get('rcsb_polymer_entity_container_identifiers'):
                auth_chains = entity['rcsb_polymer_entity_container_identifiers'].get('auth_asym_ids')
                if auth_chains and len(auth_chains) > 0:
                    # Sort chains for consistent ordering and join with commas
                    chains = ','.join(sorted(auth_chains))
            
            # Get molecule description
            molecule = 'N/A'
            if entity.get('rcsb_polymer_entity'):
                molecule = entity['rcsb_polymer_entity'].get('pdbx_description', 'N/A')
            
            # Get source organism - prioritize ncbi_scientific_name
            species = 'N/A'
            if entity.get('rcsb_entity_source_organism'):
                org_list = entity['rcsb_entity_source_organism']
                if org_list and len(org_list) > 0 and org_list[0]:
                    # Try ncbi_scientific_name first (more standardized)
                    species = org_list[0].get('ncbi_scientific_name')
                    # Fallback to scientific_name if ncbi not available
                    if not species:
                        species = org_list[0].get('scientific_name', 'N/A')
            
            # If still no species, try host organism (for recombinant proteins)
            if species == 'N/A' and entity.get('rcsb_entity_host_organism'):
                host_list = entity['rcsb_entity_host_organism']
                if host_list and len(host_list) > 0 and host_list[0]:
                    species = host_list[0].get('ncbi_scientific_name', 'N/A')
            
            # Get sequence length from entity_poly
            length = 'N/A'
            mutation_count = 0
            has_mutations = 'No'
            
            if entity.get('entity_poly'):
                poly = entity['entity_poly']
                
                # Get sequence length - handle both int and None
                seq_length = poly.get('rcsb_sample_sequence_length')
                if seq_length is not None:
                    length = seq_length
                
                # Get mutation count (this is pre-computed by RCSB)
                mut_count = poly.get('rcsb_mutation_count')
                if mut_count is not None and mut_count > 0:
                    mutation_count = mut_count
                    has_mutations = 'Yes'
            
            entities.append({
                'chains': chains,
                'molecule': molecule,
                'species': species,
                'length': length,
                'mutations': has_mutations,
                'mutation_count': mutation_count
            })
        
        return entities
        
    except requests.exceptions.RequestException as e:
        print(f"    Error fetching polymer entities for {pdb_id}: {e}")
        return []

def get_complete_pdb_info(pdb_id: str) -> List[Dict]:
    """
    Combine metadata and polymer entity information for a PDB entry.
    
    Returns one row per polymer entity, as a single PDB file can contain
    multiple proteins from different organisms.
    
    Args:
        pdb_id (str): PDB identifier
        
    Returns:
        list: List of dictionaries, one per polymer entity
    """
    # Get basic metadata
    metadata = get_pdb_metadata(pdb_id)
    if not metadata:
        return []
    
    # Small delay between API calls to be respectful
    sleep(0.2)
    
    # Get polymer entity information
    entities = get_polymer_entities(pdb_id)
    
    # If no entities found, return basic metadata only
    if not entities:
        return [{
            **metadata,
            'chains': 'N/A',
            'molecule': 'N/A',
            'species': 'N/A',
            'length': 'N/A',
            'mutations': 'N/A',
            'mutation_count': 0
        }]
    
    # Combine metadata with each entity
    results = []
    for entity in entities:
        results.append({
            **metadata,
            **entity
        })
    
    return results

def process_smiles_list(smiles_list: List[str]) -> List[Tuple]:
    """
    Process a list of SMILES strings and search PDB for each with metadata.
    
    Args:
        smiles_list (list): List of SMILES strings
        
    Returns:
        list: List of tuples containing (SMILES, metadata_dict) pairs
    """
    results = []
    
    for i, smiles in enumerate(smiles_list, 1):
        smiles = smiles.strip()
        if smiles:
            print(f"\nProcessing SMILES {i}/{len(smiles_list)}: {smiles}")
            sleep(1)  # Rate limiting for search API
            
            pdb_ids = search_pdb_by_smiles(smiles)
            if pdb_ids:
                print(f"Found {len(pdb_ids)} exact matches")
                
                for j, pdb_id in enumerate(pdb_ids, 1):
                    print(f"  Fetching metadata for {pdb_id} ({j}/{len(pdb_ids)})...")
                    
                    pdb_info_list = get_complete_pdb_info(pdb_id)
                    
                    for pdb_info in pdb_info_list:
                        results.append((smiles, pdb_info))
                    
                    # Rate limiting between PDB entries
                    if j < len(pdb_ids):
                        sleep(0.5)
            else:
                print("No exact matches found")
    
    return results

def parse_pasted_input(input_text: str) -> List[str]:
    """
    Parse input text that might come from Excel copy-paste.
    Handles both tab-separated and newline-separated input.
    
    Args:
        input_text (str): The pasted input text
        
    Returns:
        list: List of SMILES strings
    """
    lines = input_text.strip().split('\n')
    
    smiles_list = []
    for line in lines:
        # Split by tabs and take first column
        smiles = line.split('\t')[0].strip()
        if smiles:
            smiles_list.append(smiles)
    
    return smiles_list

def main():
    print("Paste your SMILES strings (can be copied directly from Excel).")
    print("Press Enter, then Ctrl+D (Unix) or Ctrl+Z (Windows) when finished:")
    
    try:
        input_text = ''
        while True:
            try:
                line = input()
                input_text += line + '\n'
            except EOFError:
                break
    except KeyboardInterrupt:
        print("\nInput cancelled.")
        return
    
    smiles_list = parse_pasted_input(input_text)
    
    if not smiles_list:
        print("No valid SMILES strings provided.")
        return
    
    print(f"\nFound {len(smiles_list)} SMILES strings to process...")
    
    # Process the list
    results = process_smiles_list(smiles_list)
    
    # Save results to CSV file
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    filename = f'pdb_matches_{timestamp}.csv'
    
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow([
            'SMILES', 
            'PDB', 
            'Chains',
            'Molecule', 
            'Species', 
            'Method', 
            'Resolution', 
            'Mutations', 
            'Mutation_Count', 
            'Length', 
            'DB_link', 
            'Title'
        ])
        
        # Write data
        for smiles, info in results:
            writer.writerow([
                smiles,
                info.get('pdb_id', 'N/A'),
                info.get('chains', 'N/A'),
                info.get('molecule', 'N/A'),
                info.get('species', 'N/A'),
                info.get('method', 'N/A'),
                info.get('resolution', 'N/A'),
                info.get('mutations', 'N/A'),
                info.get('mutation_count', 0),
                info.get('length', 'N/A'),
                info.get('db_link', 'N/A'),
                info.get('title', 'N/A')
            ])
    
    # Print summary
    print("\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)
    
    current_smiles = None
    pdb_entries = []
    unique_pdb_ids = set()
    
    for smiles, info in results:
        if smiles != current_smiles:
            if current_smiles is not None:
                print(f"\nFound {len(unique_pdb_ids)} unique PDB entries")
                print("-"*80)
            
            current_smiles = smiles
            pdb_entries = [info]
            unique_pdb_ids = {info['pdb_id']}
            print(f"\nSMILES: {smiles}")
        else:
            pdb_entries.append(info)
            unique_pdb_ids.add(info['pdb_id'])
        
        # Print abbreviated info for each entity (including chains)
        chain_info = f" [chains: {info['chains']}]" if info['chains'] != 'N/A' else ""
        print(f"  {info['pdb_id']}{chain_info}: {info['molecule'][:50]}... ({info['species']})")
    
    if current_smiles is not None:
        print(f"\nFound {len(unique_pdb_ids)} unique PDB entries")
    
    print("\n" + "="*80)
    print(f"Results saved to {filename}")
    print(f"Total rows (including multiple chains): {len(results)}")
    print("="*80)

if __name__ == "__main__":
    main()