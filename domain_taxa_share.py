#!/usr/bin/env python3

import os
import sys
import csv
import re
from collections import defaultdict
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def extract_domain_number(domain_name):
    match = re.search(r'(\d+)', domain_name)
    return int(match.group(1)) if match else float('inf')

def process_taxonomy_files(input_directory):
    domain_taxonomy = defaultdict(lambda: defaultdict(int))
    domain_totals = defaultdict(int)

    for filename in os.listdir(input_directory):
        if filename.endswith(".tsv"):
            domain = filename.split('.')[0]
            with open(os.path.join(input_directory, filename), 'r', newline='') as tsvfile:
                reader = csv.reader(tsvfile, delimiter='\t')
                next(reader)
                for row in reader:
                    if len(row) == 2:
                        cluster_id, taxonomy = row
                        domain_taxonomy[domain][taxonomy] += 1
                        domain_totals[domain] += 1

    return domain_taxonomy, domain_totals

def calculate_percentages(domain_taxonomy, domain_totals):
    domain_percentages = defaultdict(dict)

    for domain, taxonomies in domain_taxonomy.items():
        total = domain_totals[domain]
        for taxonomy, count in taxonomies.items():
            percentage = (count / total) * 100
            domain_percentages[domain][taxonomy] = percentage

    return domain_percentages

def create_stacked_chart(domain_percentages, output_svg, threshold=0.5):
    df = pd.DataFrame(domain_percentages).T
    df = df.sort_index(key=lambda x: x.map(extract_domain_number))

    # Filter out taxa that don't pass the threshold in any column
    df = df.loc[:, (df > threshold).any()]

    # Sort columns by total percentage across all domains
    column_sums = df.sum().sort_values(ascending=False)
    df = df[column_sums.index]

    # Create a color map for all taxa
    n_colors = len(df.columns)
    colors = plt.cm.get_cmap('tab20')(np.linspace(0, 1, n_colors))

    ax = df.plot(kind='bar', stacked=True, figsize=(5.5, 7.4), color=colors)

    plt.title('', fontsize=16)
    plt.xlabel('Domains', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.legend(title='Taxa', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    plt.savefig(output_svg, format='svg', dpi=300, bbox_inches='tight')
    
    print(f"Stacked chart saved as {output_svg}")

def main():
    parser = argparse.ArgumentParser(description="Generate taxonomic distribution stacked chart from taxonomy files.")
    parser.add_argument("input_directory", help="Directory containing taxonomy TSV files")
    parser.add_argument("output_svg", help="Output SVG file for the stacked chart")
    parser.add_argument("--threshold", type=float, default=0.5, help="Threshold percentage for displaying taxa (default: 0.5)")
    
    args = parser.parse_args()

    if not os.path.isdir(args.input_directory):
        print(f"Error: Input directory '{args.input_directory}' does not exist.")
        sys.exit(1)

    domain_taxonomy, domain_totals = process_taxonomy_files(args.input_directory)
    domain_percentages = calculate_percentages(domain_taxonomy, domain_totals)
    create_stacked_chart(domain_percentages, args.output_svg, args.threshold)

if __name__ == "__main__":
    main()