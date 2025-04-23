#!/usr/bin/env python3
import sys
import errno
import re
import json
import ssl
from urllib import request
from urllib.error import HTTPError, URLError
from time import sleep
from pathlib import Path
import argparse
from datetime import datetime
import socket
import io
from tqdm import tqdm

def log_message(message, error=False):
    """Print timestamped log message."""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    output = f"[{timestamp}] {message}"
    if error:
        print(output, file=sys.stderr)
    else:
        print(output)
    sys.stdout.flush()  # Ensure immediate output

def fetch_with_timeout(url, headers, timeout=30):
    """Fetch URL with timeout and return response."""
    context = ssl._create_unverified_context()
    try:
        req = request.Request(url, headers=headers)
        return request.urlopen(req, context=context, timeout=timeout)
    except socket.timeout:
        log_message(f"Timeout after {timeout} seconds", error=True)
        raise
    except URLError as e:
        log_message(f"URL Error: {str(e)}", error=True)
        raise
    except Exception as e:
        log_message(f"Fetch error: {str(e)}", error=True)
        raise

def process_pfam(pfam_id, output_dir, pbar):
    pfam_id = pfam_id.split('.')[0]  # Remove version number
    base_url = f"https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/pfam/{pfam_id}/?page_size=100"
    
    log_message(f"Processing Pfam ID: {pfam_id}")
    output_file = Path(output_dir) / f"{pfam_id}_taxonomy.txt"
    
    records_processed = 0
    page_count = 0
    next_url = base_url
    
    # Open file with line buffering
    with open(output_file, 'w', buffering=1) as f:
        f.write("protein_accession\ttax_id\n")
        
        with tqdm(desc=f"Processing {pfam_id}", unit=" records", leave=False) as pfam_pbar:
            while next_url:
                attempts = 0
                while attempts < 3:
                    try:
                        log_message(f"Fetching page {page_count + 1}")
                        
                        headers = {
                            "Accept": "application/json",
                            "User-Agent": "Python/3.x Data Collection Script"
                        }
                        
                        response = fetch_with_timeout(next_url, headers)
                        content = response.read()
                        log_message(f"Received {len(content)} bytes")
                        
                        payload = json.loads(content)
                        next_url = payload.get("next")
                        
                        # Process records and force flush after each batch
                        results = payload.get("results", [])
                        for item in results:
                            metadata = item["metadata"]
                            f.write(f"{metadata['accession']}\t{metadata['source_organism']['taxId']}\n")
                        f.flush()  # Force flush after each batch
                        
                        records_processed += len(results)
                        page_count += 1
                        log_message(f"Processed {len(results)} records on page {page_count}")
                        
                        pfam_pbar.update(len(results))
                        pbar.update(len(results))
                        
                        break
                        
                    except (socket.timeout, URLError) as e:
                        attempts += 1
                        log_message(f"Attempt {attempts}/3 failed: {str(e)}", error=True)
                        if attempts < 3:
                            sleep(10)
                        else:
                            log_message(f"Failed to process page after 3 attempts", error=True)
                            raise
                    except Exception as e:
                        log_message(f"Unexpected error: {str(e)}", error=True)
                        raise
                
                if next_url:
                    sleep(1)
    
    log_message(f"Completed {pfam_id}: processed {records_processed} records across {page_count} pages")
    return records_processed

def main():
    parser = argparse.ArgumentParser(description='Extract protein accessions and taxonomy IDs for Pfam entries')
    parser.add_argument('input_file', help='File containing list of Pfam IDs')
    parser.add_argument('output_dir', help='Directory for output files')
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(args.input_file, 'r') as f:
        pfam_ids = [line.strip() for line in f if line.strip()]
    
    log_message(f"Found {len(pfam_ids)} Pfam IDs to process")
    
    total_records = 0
    successful_pfams = 0
    failed_pfams = []
    
    with tqdm(total=len(pfam_ids), desc="Overall Progress", unit=" Pfam IDs") as pbar:
        for i, pfam_id in enumerate(pfam_ids, 1):
            log_message(f"Starting Pfam ID {i}/{len(pfam_ids)}: {pfam_id}")
            try:
                records = process_pfam(pfam_id, output_dir, pbar)
                total_records += records
                successful_pfams += 1
            except Exception as e:
                log_message(f"Failed to process {pfam_id}: {str(e)}", error=True)
                failed_pfams.append(pfam_id)
                continue
            finally:
                pbar.update(1)
    
    log_message("\nProcessing Summary:")
    log_message(f"Successfully processed: {successful_pfams}/{len(pfam_ids)} Pfam IDs")
    log_message(f"Total records processed: {total_records}")
    if failed_pfams:
        log_message("Failed Pfam IDs:", error=True)
        for pfam in failed_pfams:
            log_message(f"  - {pfam}", error=True)

if __name__ == "__main__":
    main()
