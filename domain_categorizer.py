#!/usr/bin/env python3
import argparse
import numpy as np
import sys
import re

def parse_protein_lengths(file_path):
    """
    Parse a tab-delimited file with protein IDs and lengths.
    Returns a dictionary mapping protein IDs to their lengths.
    """
    protein_lengths = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                try:
                    protein_id = parts[0]
                    length = int(parts[1])
                    protein_lengths[protein_id] = length
                except ValueError:
                    print(f"Warning: Could not parse length for protein: {line}", file=sys.stderr)
    
    return protein_lengths

def parse_data(file_path):
    clusters = []
    cluster = None
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cluster is not None:
                    clusters.append(cluster)
                try:
                    # Extract cluster size from header
                    match = re.search(r'\[(\d+)\]$', line)
                    if match:
                        cluster_size = int(match.group(1))
                        header_no_arrow = line.lstrip(">").split("[")[0]
                        representative = header_no_arrow.split("_")[0]
                        cluster = {
                            'representative': representative,
                            'size': cluster_size,
                            'members': [],
                            'original_header': line  # Store the original header
                        }
                    else:
                        print(f"Warning: Could not extract cluster size from header: {line}", file=sys.stderr)
                        cluster = None
                except Exception as e:
                    print(f"Error processing header line: {line}\n{e}", file=sys.stderr)
                    cluster = None
            elif cluster is not None:
                cluster['members'].append(line)
    
    if cluster is not None:
        clusters.append(cluster)
    return clusters

def calculate_cluster_metrics(clusters, protein_lengths, max_coverage=0.6, max_filtered_ratio=0.3, domain_size_limit=600):
    metrics = []
    # Regex to match the coordinates at the end of gene IDs
    coord_pattern = re.compile(r'(.+)_(\d+)_(\d+)(?:\s*$|\[)')
    
    # Set to track missing protein IDs
    missing_ids = set()
    total_filtered = 0
    total_members = 0
    clusters_filtered = 0
    
    for cluster in clusters:
        domain_sizes = []
        avg_coords = []
        protein_sizes = []  # Store protein sizes for this cluster
        coords_with_sizes = []  # Store coordinates only when we have a protein size
        filtered_count = 0
        valid_members = []  # Store valid members after filtering
        total_members += len(cluster['members'])
        
        for member in cluster['members']:
            match = coord_pattern.search(member)
            if not match:
                print(f"Skipping member with unexpected format: {member}", file=sys.stderr)
                continue
            
            try:
                protein_base = match.group(1)  # Everything before the last _NUM_NUM
                coord1 = int(match.group(2))
                coord2 = int(match.group(3))
                
                # Ensure coord2 is greater than coord1
                if coord2 < coord1:
                    coord1, coord2 = coord2, coord1
                
                domain_size = coord2 - coord1
                average_coord = (coord1 + coord2) / 2.0
                
                # Try different variations of the protein ID to find a match
                protein_id = None
                protein_length = None
                
                # 1. Try exact match
                if protein_base in protein_lengths:
                    protein_id = protein_base
                    protein_length = protein_lengths[protein_id]
                
                # 2. Try removing last underscore part if present
                elif '_' in protein_base:
                    parts = protein_base.split('_')
                    potential_id = '_'.join(parts[:-1])
                    if potential_id in protein_lengths:
                        protein_id = potential_id
                        protein_length = protein_lengths[potential_id]
                
                # 3. Try just the first part (before any underscore)
                elif protein_base.split('_')[0] in protein_lengths:
                    protein_id = protein_base.split('_')[0]
                    protein_length = protein_lengths[protein_id]
                
                # If we found a protein length, check coverage
                if protein_length is not None:
                    coverage = domain_size / protein_length
                    
                    # If coverage exceeds threshold, skip this member
                    if coverage > max_coverage:
                        print(f"Filtering out member with high coverage ({coverage:.2f}): {member}", file=sys.stderr)
                        filtered_count += 1
                        continue
                    
                    # Store the protein length and corresponding coordinate
                    protein_sizes.append(protein_length)
                    coords_with_sizes.append(average_coord)
                else:
                    # If protein length is not available, use the domain size limit
                    if domain_size > domain_size_limit:
                        print(f"Skipping unusually large domain size ({domain_size}) for member with unknown length: {member}", file=sys.stderr)
                        filtered_count += 1
                        continue
                    
                    # Add to missing IDs set
                    missing_ids.add(protein_base)
                
                domain_sizes.append(domain_size)
                avg_coords.append(average_coord)
                valid_members.append(member)  # Store valid member
                
            except ValueError as e:
                print(f"Error converting coordinates for {member}: {e}", file=sys.stderr)
                continue
        
        total_filtered += filtered_count
        
        # Skip clusters where too many members were filtered out
        if len(cluster['members']) > 0 and filtered_count / len(cluster['members']) > max_filtered_ratio:
            print(f"Skipping cluster with representative {cluster['representative']} - {filtered_count}/{len(cluster['members'])} members filtered out", file=sys.stderr)
            clusters_filtered += 1
            continue
        
        if not domain_sizes or not avg_coords:
            print(f"No valid members in cluster with representative {cluster['representative']}", file=sys.stderr)
            continue
        
        avg_domain_size = np.mean(domain_sizes)
        overall_avg_coord = np.mean(avg_coords)
        
        # Calculate average protein size and relative position for this cluster
        if protein_sizes and coords_with_sizes:
            avg_protein_size = np.mean(protein_sizes)
            avg_coord_with_size = np.mean(coords_with_sizes)
            
            # Calculate relative position (normalized by protein size)
            # Clamp the value between 0 and 1
            relative_position = min(1.0, max(0.0, avg_coord_with_size / avg_protein_size))
        else:
            # If no protein sizes available, use a default value
            print(f"Warning: No protein lengths available for cluster with representative {cluster['representative']}, using default relative position", file=sys.stderr)
            avg_protein_size = 0
            relative_position = 0.5  # Default to middle of color range
        
        metrics.append({
            'cluster_size': cluster['size'],
            'avg_domain_size': avg_domain_size,
            'avg_coordinate': overall_avg_coord,
            'avg_protein_size': avg_protein_size,
            'relative_position': relative_position,
            'original_header': cluster['original_header'],
            'valid_members': valid_members  # Store only the valid members after filtering
        })
    
    # Report filtering statistics
    print(f"Filtering summary: {total_filtered}/{total_members} members ({total_filtered/total_members:.1%}) and {clusters_filtered}/{len(clusters)} clusters ({clusters_filtered/len(clusters):.1%}) filtered out", file=sys.stderr)
    
    # Report missing protein IDs
    if missing_ids:
        print(f"Warning: {len(missing_ids)} protein IDs were not found in the protein lengths file.", file=sys.stderr)
        if len(missing_ids) <= 10:
            print("Missing IDs:", list(missing_ids), file=sys.stderr)
        else:
            print("First 10 missing IDs:", list(missing_ids)[:10], file=sys.stderr)
    
    return metrics

def write_domain_files(metrics, output_prefix):
    """
    Write three output files for N-terminal, middle, and C-terminal domains.
    Categorize based on relative position:
    - N-terminal: 0-0.3
    - Middle: 0.3-0.7
    - C-terminal: 0.7-1.0
    """
    # Prepare the output files
    n_terminal_file = f"{output_prefix}_n_terminal.txt"
    middle_file = f"{output_prefix}_middle.txt"
    c_terminal_file = f"{output_prefix}_c_terminal.txt"
    
    # Categorize the clusters
    n_terminal_clusters = []
    middle_clusters = []
    c_terminal_clusters = []
    
    for metric in metrics:
        rel_pos = metric['relative_position']
        
        if rel_pos < 0.3:
            n_terminal_clusters.append(metric)
        elif rel_pos < 0.7:
            middle_clusters.append(metric)
        else:
            c_terminal_clusters.append(metric)
    
    # Write the files
    def write_file(filename, clusters):
        with open(filename, 'w') as f:
            for cluster in clusters:
                f.write(f"{cluster['original_header']}\n")
                for member in cluster['valid_members']:
                    f.write(f"{member}\n")
                f.write("\n")  # Add empty line between clusters
    
    write_file(n_terminal_file, n_terminal_clusters)
    write_file(middle_file, middle_clusters)
    write_file(c_terminal_file, c_terminal_clusters)
    
    # Report results
    print(f"N-terminal domains ({len(n_terminal_clusters)} clusters): {n_terminal_file}")
    print(f"Middle domains ({len(middle_clusters)} clusters): {middle_file}")
    print(f"C-terminal domains ({len(c_terminal_clusters)} clusters): {c_terminal_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Categorize protein domains into N-terminal, middle, and C-terminal regions."
    )
    parser.add_argument("input_file", help="Path to the input data file")
    parser.add_argument("--output", help="Output file prefix for the three domain files", default="domains")
    parser.add_argument("--lengths", help="Path to a TSV file with protein IDs and lengths", required=True)
    parser.add_argument("--max-coverage", type=float, default=0.6, 
                        help="Maximum allowed domain coverage (domain size / protein length). Default: 0.6")
    parser.add_argument("--max-filtered-ratio", type=float, default=0.3, 
                        help="Maximum allowed ratio of filtered members in a cluster. Default: 0.3")
    parser.add_argument("--domain-size-limit", type=int, default=600,
                        help="Maximum allowed domain size for proteins without length information. Default: 600")
    args = parser.parse_args()
    
    protein_lengths = {}
    print(f"Parsing protein lengths from {args.lengths}...", file=sys.stderr)
    protein_lengths = parse_protein_lengths(args.lengths)
    print(f"Loaded lengths for {len(protein_lengths)} proteins", file=sys.stderr)
    
    print(f"Parsing data from {args.input_file}...", file=sys.stderr)
    clusters = parse_data(args.input_file)
    
    if not clusters:
        print("No cluster data found. Exiting.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(clusters)} clusters. Calculating metrics with filtering...", file=sys.stderr)
    metrics = calculate_cluster_metrics(clusters, protein_lengths, 
                                       max_coverage=args.max_coverage, 
                                       max_filtered_ratio=args.max_filtered_ratio,
                                       domain_size_limit=args.domain_size_limit)
    
    if not metrics:
        print("No valid cluster metrics computed after filtering. Exiting.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Writing {len(metrics)} clusters after filtering to domain files...", file=sys.stderr)
    write_domain_files(metrics, args.output)

if __name__ == "__main__":
    main()

