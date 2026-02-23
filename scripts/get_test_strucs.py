"""
get_test_structures.py
==============================
Fetches orthologs of a reference protein via OrthoDB, maps them to UniProt
accessions using direct cross-reference lookups, then retrieves AlphaFold 
structures where available.

Outputs
-------
  orthodb_orthologs.csv    -- OrthoDB-level data for every ortholog resolved
  alphafold_metadata.csv   -- AlphaFold metadata for structures that were found
  af_structures/           -- Downloaded PDB files (CIF converted via obabel)
"""

import csv
import subprocess
import sys
import time
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
OUTPUT_DIR      = Path(".")
AF_STRUCT_DIR   = OUTPUT_DIR / "af_structures"
ORTHODB_CSV     = OUTPUT_DIR / "orthodb_orthologs.csv"
ALPHAFOLD_CSV   = OUTPUT_DIR / "alphafold_metadata.csv"

ORTHODB_BASE    = "https://data.orthodb.org/v12"
ALPHAFOLD_BASE  = "https://alphafold.ebi.ac.uk/api"
UNIPROT_BASE    = "https://rest.uniprot.org/uniprotkb"

RATE_LIMIT_DELAY   = 0.4    # seconds between API calls
MAX_ORTHOLOGS      = 500    # safety cap -- orthogroups can be very large


def get_json(url, params=None, retries=3):
    """
    GET a JSON endpoint with exponential-backoff retry.
    Returns parsed JSON on success, None on all failures.
    404 returns None silently (expected for missing entries).
    """
    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(url, params=params, timeout=30)
            if resp.status_code == 404:
                return None
            resp.raise_for_status()
            return resp.json()
        except requests.RequestException as exc:
            wait = 2 ** attempt
            print(f"  [WARN] Attempt {attempt}/{retries} -- {url}: {exc}")
            if attempt < retries:
                time.sleep(wait)
    return None


def prompt_for_inputs():
    """Collect and validate a UniProt accession and optional taxon level."""
    print("\n============================================================")
    print("  Ortholog Discovery + AlphaFold Structure Retrieval")
    print("============================================================")
    cwd = input(f"\nPlease enter an output directory for downloaded structures:")
    if cwd.strip():
        global OUTPUT_DIR, AF_STRUCT_DIR, ORTHODB_CSV, ALPHAFOLD_CSV
        OUTPUT_DIR = Path(cwd.strip())
        AF_STRUCT_DIR = OUTPUT_DIR / "af_structures"
        ORTHODB_CSV = OUTPUT_DIR / "orthodb_orthologs.csv"
        ALPHAFOLD_CSV = OUTPUT_DIR / "alphafold_metadata.csv"
        print(f"  Outputs will be saved to: {OUTPUT_DIR.resolve()}")
    else:
        print(f"  Using current directory for outputs: {OUTPUT_DIR.resolve()}")
    
    uid = ""
    while not uid:
        uid = input(
            "\nEnter the UniProt ID of the reference protein (e.g. P10275): "
        ).strip().upper()
        if not uid:
            print("  Please enter a non-empty accession.")

    # OrthoDB uses NCBI taxonomy levels: common options include
    # 2759 (Eukaryota), 33208 (Metazoa), 7742 (Vertebrata), 40674 (Mammalia)
    print(
        "\nTaxonomic scope for ortholog search (NCBI taxon ID)."
        "\n  Examples: 2759=Eukaryota  33208=Metazoa  7742=Vertebrata"
        "\n  Leave blank to use OrthoDB's default (all taxa)."
    )
    taxon_input = input("  Taxon ID [blank = all]: ").strip()
    taxon_id = taxon_input if taxon_input.isdigit() else None

    print()
    protein_code = ""
    while not protein_code:
        protein_code = input(
            "\nEnter a short protein code for output filenames (e.g. AR, BRCA1): "
        ).strip()
        if not protein_code:
            print("  Please enter a non-empty protein code.")

    print()
    return uid, taxon_id, protein_code


def genesearch(uniprot_id):
    """
    /genesearch returns a single best-matching gene entry.
    We extract the OrthoDB gene param (e.g. '9606_1:0052ae') which is the
    handle needed for all subsequent /gene and /orthologs calls.

    Response shape (relevant fields):
      {"gene": {"gene_id": {"param": "9606_1:0052ae"}}}
    """
    print(f"  Querying /genesearch for UniProt ID: {uniprot_id}")
    data = get_json(f"{ORTHODB_BASE}/genesearch", params={"query": uniprot_id})
    if not data:
        print(f"  [ERROR] /genesearch returned no data for {uniprot_id}.")
        return None
    gene = data.get("gene", {})
    gene_id_obj = gene.get("gene_id", {})
    gene_param = gene_id_obj.get("param")
    if not gene_param:
        print(f"  [ERROR] Could not find gene_id.param for {uniprot_id}.")
        return None
    print(f"  Found OrthoDB gene param: {gene_param}")
    return gene_param


def get_orthologs(gene_param, taxon_id=None):
    """
    /orthologs returns one record per ortholog gene.

    Response shape (relevant fields):
      {"data": [{"gene": {"id": "...", "param": "392033_0:007e62"},
                 "taxon_id": "392033_0",
                 "clade_id": 33208}, ...]}

    taxon_id values from OrthoDB are formatted as "392033_0" where the numeric
    prefix before the underscore is the NCBI taxonomy ID. We split on "_" and
    take the first element to obtain a clean integer taxon ID for downstream
    NCBI taxonomy lookups.

    We store both the OrthoDB gene param (for the /gene xref lookup in Step 4)
    and available taxonomy data for the CSV.  The 'id' field is a readable name
    while 'param' is the stable API handle.
    """
    print(f"  Querying /orthologs for gene param: {gene_param}")
    params = {"id": gene_param}
    if taxon_id:
        params["taxon"] = taxon_id
    data = get_json(f"{ORTHODB_BASE}/orthologs", params=params)
    if not data or "data" not in data:
        print(f"  [ERROR] No ortholog data returned for {gene_param}.")
        return []

    orthologs = []
    for entry in data["data"]:
        gene_obj = entry.get("gene", {})
        param = gene_obj.get("param")
        readable_id = gene_obj.get("id", "")
        raw_taxon = entry.get("taxon_id", "")
        clade = entry.get("clade_id", "")

        # Split "392033_0" on "_" and take the first element to get the
        # bare NCBI taxon ID integer string (e.g. "392033").
        ncbi_taxon_id = raw_taxon.split("_")[0] if raw_taxon else ""

        if param:
            orthologs.append({
                "orthodb_param":  param,
                "orthodb_id":     readable_id,
                "taxon_id":       ncbi_taxon_id,   # clean NCBI taxon ID
                "raw_taxon_id":   raw_taxon,        # preserved for reference
                "clade_id":       clade,
            })

    # Deduplicate on orthodb_param; a gene can appear in multiple clades
    seen = set()
    unique = []
    for o in orthologs:
        if o["orthodb_param"] not in seen:
            seen.add(o["orthodb_param"])
            unique.append(o)

    total = len(data["data"])
    print(
        f"  Found {total} ortholog entries, {len(unique)} unique after dedup."
    )

    return unique


def get_uniprot(orthodb_param, orthodb_id):
    query_term = orthodb_id if orthodb_id else orthodb_param
    data = get_json(
        f"{UNIPROT_BASE}/search",
        params={"query": query_term, "format": "json", "size": 1},
    )
    if not data or not data.get("results"):
        return None
    hit = data["results"][0]
    organism = hit.get("organism", {})
    return {
        "accession":       hit.get("primaryAccession"),
        "scientific_name": organism.get("scientificName", ""),
        "taxon_id":        str(organism.get("taxonId", "")),
    }


def resolve_uniprot_id(ortholog):
    param = ortholog["orthodb_param"]
    hit = get_uniprot(param, ortholog.get("orthodb_id", ""))
    if hit and hit["accession"]:
        if hit["taxon_id"] == str(ortholog["taxon_id"]):
            ortholog["scientific_name"]    = hit["scientific_name"]
            ortholog["search_result"]  = "success"
            return hit["accession"]
        else:
            print(f"    [MISMATCH] {param}: UniProt taxon {hit['taxon_id']} "
                  f"!= OrthoDB taxon {ortholog['taxon_id']} — rejecting.")
    ortholog["search_result"] = "failed"
    return None


def fetch_alphafold_entries(uniprot_id):
    """
    The /uniprot/summary/{id}.json endpoint returns all structures associated
    with a UniProt accession, which can include multiple segments for long
    proteins.

    We iterate all entries rather than only structures[0] so that multi-segment
    proteins (>2700 aa) get all their fragments captured.

    Response shape (relevant fields):
      {
        "uniprot_entry": {"sequence_length": 919},
        "structures": [
          {"summary": {
              "model_identifier": "AF-P10275-F1-model_v4",
              "model_url": "https://...",
              "model_category": "AB-INITIO",
              "provider": "AlphaFold DB",
              "confidence_type": "pLDDT",
              "confidence_avg_local_score": 84.2,
              "entities": [{"identifier": "P10275",
                             "description": "Androgen receptor", ...}]
          }}
        ]
      }
    """
    data = get_json(f"{ALPHAFOLD_BASE}/uniprot/summary/{uniprot_id}.json")
    if not data or not data.get("structures"):
        return []

    seq_length = data.get("uniprot_entry", {}).get("sequence_length", "")
    results = []

    for struct in data["structures"]:
        s = struct.get("summary", {})
        model_id  = s.get("model_identifier", "")
        model_url = s.get("model_url", "")
        if not model_id or not model_url:
            continue

        entity = s.get("entities", [{}])[0]
        results.append({
            "model_id":          model_id,
            "model_url":         model_url,
            "model_category":    s.get("model_category", ""),
            "provider":          s.get("provider", ""),
            "confidence_type":   s.get("confidence_type", ""),
            "confidence_score":  s.get("confidence_avg_local_score", ""),
            "sequence_length":   seq_length,
            "description":       entity.get("description", ""),
        })

    return results


def download_structure(model_url, output_path):
    """Download a structure file to output_path. Returns True on success."""
    try:
        resp = requests.get(model_url, timeout=60)
        resp.raise_for_status()
        output_path.write_bytes(resp.content)
        return True
    except requests.RequestException as exc:
        print(f"  [ERROR] Download failed for {model_url}: {exc}")
        return False


def convert_cif_to_pdb(cif_path):
    """
    Convert a CIF file to PDB format using Open Babel (obabel).

    Open Babel is invoked as a subprocess:
        obabel input.cif -O output.pdb

    The -O flag specifies the output file; obabel infers the output format
    from the extension. We do not pass -h (add hydrogens) or any other
    modification flags because we want to preserve the original coordinates
    exactly as deposited.

    On success:
      - The new .pdb file exists on disk
      - The original .cif file is deleted to avoid accumulating redundant files
      - Returns the Path to the new .pdb file

    On failure (obabel not found, non-zero exit code, or output not created):
      - Prints a warning and returns the original cif_path unchanged so the
        caller can still record the CIF file in the metadata CSV rather than
        losing track of it entirely.
    """
    pdb_path = cif_path.with_suffix(".pdb")
    try:
        result = subprocess.run(
            ["obabel", str(cif_path), "-O", str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0 or not pdb_path.exists():
            print(f"    [WARN] obabel conversion failed for {cif_path.name}.")
            print(f"           stderr: {result.stderr.strip()}")
            return cif_path   # fall back to keeping CIF
        # Conversion succeeded -- clean up the intermediate CIF
        cif_path.unlink()
        print(f"    Converted CIF → PDB: {pdb_path.name}")
        return pdb_path
    except FileNotFoundError:
        print(
            "    [WARN] obabel not found on PATH. "
            "CIF file retained without conversion."
        )
        return cif_path
    except subprocess.TimeoutExpired:
        print(f"    [WARN] obabel timed out for {cif_path.name}. CIF retained.")
        return cif_path


# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

ORTHODB_FIELDS = [
    "orthodb_param", "orthodb_id", "taxon_id", "scientific_name",
    "clade_id", "uniprot_id", "search_result",
]

ALPHAFOLD_FIELDS = [
    "uniprot_id", "orthodb_param", "taxon_id", "scientific_name",
    "model_id", "model_url", "model_category", "provider",
    "confidence_type", "confidence_score",
    "sequence_length", "description", "structure_file",
]


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved {len(rows)} rows → {path}")


def main():
    """
    Main execution flow:
    1. Prompt user for UniProt ID, taxon ID, and protein code.
    2. Use /genesearch to find the OrthoDB gene param for the reference protein.
    3. Use /orthologs to retrieve orthologs, filtering by taxon if specified.
    3b. For each unique taxon ID among the orthologs, query NCBI taxonomy to get
        the scientific name and populate it on the ortholog records.
    4. For each ortholog, attempt to resolve a UniProt ID via the /gene xref
        lookup. Store the resolution method and result on the ortholog record.
    5. For each resolved UniProt ID, query the AlphaFold API for associated
        structures. Download each structure file, convert from CIF to PDB if
        needed, and rename to a consistent "[Genus-species]_[protein_code].pdb" format.
        Save OrthoDB-level data for all orthologs to orthodb_orthologs.csv and
        AlphaFold metadata for resolved structures to alphafold_metadata.csv.
    """
    # ------------------------------------------------------------------
    # Step 1: User input
    # ------------------------------------------------------------------
    reference_uid, taxon_id, protein_code = prompt_for_inputs()

    # ------------------------------------------------------------------
    # Step 2: Locate the reference gene in OrthoDB
    # ------------------------------------------------------------------
    gene_param = genesearch(reference_uid)
    if not gene_param:
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 3: Retrieve orthologs
    # ------------------------------------------------------------------
    orthologs = get_orthologs(gene_param, taxon_id)
    if not orthologs:
        print("  No orthologs found. Exiting.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # Step 3b: Set ortholog limit
    # ------------------------------------------------------------------
    global MAX_ORTHOLOGS
    print(f"\n  {len(orthologs)} ortholog(s) found.")
    set_max = input(
        f"  Enter a number to limit orthologs, 'all' to use all {len(orthologs)}, "
        f"or press Enter to use the default ({MAX_ORTHOLOGS}): "
    ).strip()
    if set_max.isdigit():
        MAX_ORTHOLOGS = int(set_max)
    elif set_max == "all":
        MAX_ORTHOLOGS = len(orthologs)
    else:
        print(f"  Keeping default limit of {MAX_ORTHOLOGS} orthologs.")
    if len(orthologs) > MAX_ORTHOLOGS:
        print(f"  Limiting to {MAX_ORTHOLOGS} orthologs.")
        orthologs = orthologs[:MAX_ORTHOLOGS]

    # ------------------------------------------------------------------
    # Step 4: Resolve each ortholog to a UniProt accession
    # ------------------------------------------------------------------
    print(f"\n[Step 4] Resolving {len(orthologs)} ortholog(s) to UniProt IDs...")
    for i, orth in enumerate(orthologs, 1):
        uid = resolve_uniprot_id(orth)
        orth["uniprot_id"] = uid or ""
        status = uid if uid else "-- not resolved --"
        print(f"  [{i}/{len(orthologs)}] {orth['orthodb_param']} → {status}"
              f"  (method: {orth.get('search_result', 'n/a')})")
        time.sleep(RATE_LIMIT_DELAY)

    # Save OrthoDB-level data regardless of AlphaFold availability.
    # scientific_name is included via ORTHODB_FIELDS.
    write_csv(ORTHODB_CSV, ORTHODB_FIELDS, orthologs)

    resolved = [o for o in orthologs if o["uniprot_id"]]
    print(f"\n  Resolved {len(resolved)}/{len(orthologs)} orthologs to UniProt IDs.")

    # ------------------------------------------------------------------
    # Step 5: Fetch AlphaFold structures and download
    # ------------------------------------------------------------------
    print(f"\n[Step 5] Querying AlphaFold DB for {len(resolved)} UniProt ID(s)...")
    AF_STRUCT_DIR.mkdir(exist_ok=True)
    af_records = []

    for i, orth in enumerate(resolved, 1):
        uid = orth["uniprot_id"]
        print(f"  [{i}/{len(resolved)}] AlphaFold lookup: {uid}")
        entries = fetch_alphafold_entries(uid)

        if not entries:
            print(f"    No AlphaFold structure found for {uid}.")
            time.sleep(RATE_LIMIT_DELAY)
            continue

        for entry in entries:
            model_id  = entry["model_id"]
            model_url = entry["model_url"]

            # Build a sanitised species slug from the scientific name.
            # "Homo sapiens" → "Homo-sapiens"; fall back to taxon_id if blank.
            sci_name = orth.get("scientific_name", "").strip()
            if sci_name:
                # Keep only genus + species (first two words), replace space
                # with hyphen, strip any characters unsafe for filenames.
                parts = sci_name.split()[:2]
                species_slug = "-".join(parts)
                species_slug = "".join(
                    c for c in species_slug if c.isalnum() or c in "-_"
                )
            else:
                species_slug = orth.get("taxon_id", uid)

            # Infer file extension from URL (usually .pdb or .cif)
            ext = Path(model_url).suffix or ".pdb"
            # Use a temporary name for download; rename after any conversion.
            tmp_path = AF_STRUCT_DIR / f"_tmp_{uid}_{model_id}{ext}"

            print(f"    Downloading {model_id} → {tmp_path.name}")
            success = download_structure(model_url, tmp_path)

            if success:
                # Convert CIF → PDB via obabel if needed.
                if tmp_path.suffix.lower() == ".cif":
                    tmp_path = convert_cif_to_pdb(tmp_path)

                # Rename to final "[Genus-species]_[protein_code].pdb".
                # If multiple fragments exist for the same protein, append
                # the fragment index (F1, F2 …) extracted from the model_id
                # to avoid collisions (e.g. "Homo-sapiens_AR_F2.pdb").
                import re as _re
                frag_match = _re.search(r"-(F\d+)-", model_id)
                frag_suffix = f"_{frag_match.group(1)}" if frag_match else ""
                final_name = f"{species_slug}_{protein_code}{frag_suffix}.pdb"
                final_path = AF_STRUCT_DIR / final_name
                # Avoid silent overwrites if two orthologs share a name
                counter = 1
                while final_path.exists():
                    final_path = AF_STRUCT_DIR / (
                        f"{species_slug}_{protein_code}{frag_suffix}"
                        f"_{counter}.pdb"
                    )
                    counter += 1
                tmp_path.rename(final_path)
                out_path = final_path
                print(f"    Renamed → {out_path.name}")

                record = {
                    "uniprot_id":      uid,
                    "orthodb_param":   orth["orthodb_param"],
                    "taxon_id":        orth["taxon_id"],
                    "scientific_name": orth["scientific_name"],
                    "structure_file":  out_path.name,
                    **entry,
                }
                af_records.append(record)
            else:
                print(f"    [WARN] Skipping {model_id} due to download failure.")

        time.sleep(RATE_LIMIT_DELAY)

    # ------------------------------------------------------------------
    # Sort AlphaFold records by confidence_score descending before saving.
    # Records where confidence_score is empty (missing) are sorted to the
    # bottom by coercing them to -1.0 so the sort key is always numeric.
    # ------------------------------------------------------------------
    if af_records:
        af_records.sort(
            key=lambda r: float(r["confidence_score"])
            if r["confidence_score"] != ""
            else -1.0,
            reverse=True,
        )
        write_csv(ALPHAFOLD_CSV, ALPHAFOLD_FIELDS, af_records)
    else:
        print("\n  No AlphaFold structures were downloaded.")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n============================================================")
    print(f"  Orthologs queried       : {len(orthologs)}")
    print(f"  UniProt IDs resolved    : {len(resolved)}")
    print(f"  AlphaFold structures    : {len(af_records)}")
    print(f"  OrthoDB CSV             : {ORTHODB_CSV}")
    if af_records:
        print(f"  AlphaFold metadata CSV  : {ALPHAFOLD_CSV}")
        print(f"  Structure files         : {AF_STRUCT_DIR}/")
    print("============================================================\n")


if __name__ == "__main__":
    main()