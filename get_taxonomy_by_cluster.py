#!/usr/bin/env python3

import os
import sys
import csv
import argparse
from collections import defaultdict

def read_tsv_file(file_path, taxonomy_column):
    encodings = ['utf-8', 'iso-8859-1', 'windows-1252']
    
    for encoding in encodings:
        try:
            with open(file_path, 'r', newline='', encoding=encoding) as tsvfile:
                reader = csv.reader(tsvfile, delimiter='\t')
                next(reader)  # Skip header row
                for row in reader:
                    if len(row) > max(6, taxonomy_column - 1):
                        yield row
            return  # If successful, exit the function
        except UnicodeDecodeError:
            continue
    
    # If all encodings fail, try binary mode
    with open(file_path, 'rb') as tsvfile:
        for line in tsvfile:
            try:
                decoded_line = line.decode('utf-8').strip()
            except UnicodeDecodeError:
                decoded_line = line.decode('iso-8859-1').strip()
            row = decoded_line.split('\t')
            if len(row) > max(6, taxonomy_column - 1):
                yield row

def process_files(input_directory, tsv_database_file, taxonomy_column):
    # Dictionary to store cluster IDs from input files
    cluster_ids = set()

    # Read cluster IDs from all input files
    for filename in os.listdir(input_directory):
        if filename.endswith("_sorted_unique"):
            with open(os.path.join(input_directory, filename), 'r') as f:
                cluster_ids.update(line.strip() for line in f)

    # Dictionary to store taxonomy information
    taxonomy_info = defaultdict(set)

    # Read TSV database file and extract taxonomy information
    for row in read_tsv_file(tsv_database_file, taxonomy_column):
        cluster_id = row[6]  # Assuming cluster ID is in the 7th column (index 6)
        if cluster_id in cluster_ids:
            tax_info = row[taxonomy_column - 1].strip()
            if tax_info:
                taxonomy_info[cluster_id].add(tax_info)

    # Write output files
    for filename in os.listdir(input_directory):
        if filename.endswith("_sorted_unique"):
            input_file = os.path.join(input_directory, filename)
            output_file = f"{input_file}_taxonomy.tsv"
            
            with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
                writer = csv.writer(outfile, delimiter='\t')
                writer.writerow(["Cluster_ID", "Taxonomic_Info"])
                
                for line in infile:
                    cluster_id = line.strip()
                    if cluster_id in taxonomy_info:
                        for tax_info in taxonomy_info[cluster_id]:
                            writer.writerow([cluster_id, tax_info])
            
            print(f"Created taxonomy file: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract taxonomy information for cluster IDs.")
    parser.add_argument("input_directory", help="Directory containing *_sorted_unique files")
    parser.add_argument("tsv_database_file", help="TSV file with taxonomic information")
    parser.add_argument("taxonomy_column", type=int, help="Column number in the TSV file containing the taxonomy information")
    
    args = parser.parse_args()

    if not os.path.isdir(args.input_directory):
        print(f"Error: Input directory '{args.input_directory}' does not exist.")
        sys.exit(1)

    if not os.path.isfile(args.tsv_database_file):
        print(f"Error: TSV database file '{args.tsv_database_file}' does not exist.")
        sys.exit(1)

    process_files(args.input_directory, args.tsv_database_file, args.taxonomy_column)

if __name__ == "__main__":
    main()