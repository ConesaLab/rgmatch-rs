#!/usr/bin/env python3
import os
import gzip
import requests
import json
import sys

# Constants
GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
ENCODE_SEARCH_URL = "https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=TF+ChIP-seq&target.label=CTCF&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assembly=GRCh38&files.file_type=bed+narrowPeak&format=json"

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "../tests/data")
GTF_OUTPUT = os.path.join(OUTPUT_DIR, "benchmark.gtf")
BED_OUTPUT = os.path.join(OUTPUT_DIR, "benchmark.bed")

VALID_CHROMOSOMES = {"chr21", "chr22"}

def ensure_output_dir():
    """Ensure the output directory exists."""
    if not os.path.exists(OUTPUT_DIR):
        print(f"Creating output directory: {OUTPUT_DIR}")
        os.makedirs(OUTPUT_DIR)

def stream_gencode_gtf():
    """Stream and filter GENCODE GTF."""
    print(f"Streaming GENCODE GTF from {GENCODE_URL}...")
    try:
        with requests.get(GENCODE_URL, stream=True) as r:
            r.raise_for_status()
            # Use gzip.open on the raw stream
            with gzip.open(r.raw, mode='rt') as f_in, open(GTF_OUTPUT, 'w') as f_out:
                print(f"Filtering for {VALID_CHROMOSOMES}...")
                count = 0
                for line in f_in:
                    if line.startswith('#'):
                        f_out.write(line) # Keep headers
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) > 0 and parts[0] in VALID_CHROMOSOMES:
                        f_out.write(line)
                        count += 1
                print(f"  -> Wrote {count} lines to {GTF_OUTPUT}")
    except Exception as e:
        print(f"Error streaming GENCODE GTF: {e}", file=sys.stderr)
        sys.exit(1)

def get_encode_bed_url():
    """Query ENCODE API for the first available bed narrowPeak file."""
    print(f"Querying ENCODE API...")
    try:
        # Step 1: Search for an Experiment
        search_url = ENCODE_SEARCH_URL + "&limit=1"
        r = requests.get(search_url, headers={"Accept": "application/json"})
        r.raise_for_status()
        data = r.json()
        
        if '@graph' not in data or len(data['@graph']) == 0:
             print("Error: No experiments found in ENCODE search.", file=sys.stderr)
             sys.exit(1)
             
        experiment = data['@graph'][0]
        experiment_id = experiment.get('@id')
        if not experiment_id:
             print("Error: Experiment ID not found.", file=sys.stderr)
             sys.exit(1)

        print(f"  -> Found experiment: {experiment_id}")
        
        # Step 2: Fetch Experiment Details to get the full file list
        experiment_url = f"https://www.encodeproject.org{experiment_id}"
        r_exp = requests.get(experiment_url, headers={"Accept": "application/json"})
        r_exp.raise_for_status()
        exp_data = r_exp.json()

        for file_obj in exp_data.get('files', []):
            if file_obj.get('file_type') == 'bed narrowPeak' and file_obj.get('status') == 'released':
                href = file_obj.get('href')
                if href:
                    full_url = f"https://www.encodeproject.org{href}"
                    print(f"  -> Found BED file: {full_url}")
                    return full_url
        
        print("Error: No suitable BED file found in this experiment.", file=sys.stderr)
        sys.exit(1)

    except Exception as e:
        print(f"Error querying ENCODE API: {e}", file=sys.stderr)
        sys.exit(1)

def stream_encode_bed(url):
    """Stream and filter ENCODE BED file."""
    print(f"Streaming ENCODE BED from {url}...")
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            # ENCODE bed files are usually gzipped
            with gzip.open(r.raw, mode='rt') as f_in, open(BED_OUTPUT, 'w') as f_out:
                print(f"Filtering for {VALID_CHROMOSOMES}...")
                count = 0
                for line in f_in:
                    parts = line.split('\t')
                    if len(parts) > 0 and parts[0] in VALID_CHROMOSOMES:
                        f_out.write(line)
                        count += 1
                
                # Edge Case Injection
                print("Injecting edge case 'Missing Chromosome'...")
                f_out.write("chr99\t100\t200\tfake_region\t0\t.\t.\t.\t.\t.\n")
                
                print(f"  -> Wrote {count} lines + 1 edge case to {BED_OUTPUT}")

    except Exception as e:
         print(f"Error streaming ENCODE BED: {e}", file=sys.stderr)
         sys.exit(1)

def main():
    ensure_output_dir()
    stream_gencode_gtf()
    bed_url = get_encode_bed_url()
    stream_encode_bed(bed_url)
    print("Done! Benchmark dataset generated.")

if __name__ == "__main__":
    main()
