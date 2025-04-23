#!/usr/bin/env python3
import argparse
import sys
from collections import Counter

def is_header_line(line):
    """Check if a line appears to be a header."""
    parts = line.strip().split()
    header_words = ['gene_oid', 'clust_ID', 'clust_ID_x']
    return any(word in line for word in header_words) and not line.strip().split()[0].isdigit()

def analyze_gene_clusters(cluster_file_path, gene_list_file_path, debug=False):
    """
    Analyze how many genes from a given list appear in each cluster.
    """
    try:
        with open(gene_list_file_path, 'r') as f:
            genes_to_check = set(line.strip() for line in f if line.strip())
        if debug:
            print(f"Loaded {len(genes_to_check)} genes to check", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: Gene list file '{gene_list_file_path}' not found!", file=sys.stderr)
        sys.exit(1)
        
    cluster_counts = {}
    
    try:
        with open(cluster_file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or is_header_line(line):
                    continue
                
                parts = line.split()
                if len(parts) < 2:
                    continue
                
                try:
                    gene_id = parts[0]
                    cluster = int(parts[1])
                    
                    if gene_id in genes_to_check:
                        cluster_counts[cluster] = cluster_counts.get(cluster, 0) + 1
                            
                except ValueError:
                    continue
                    
    except FileNotFoundError:
        print(f"Error: Cluster file '{cluster_file_path}' not found!", file=sys.stderr)
        sys.exit(1)
    
    return cluster_counts

def create_cluster_files(cluster_counts):
    """
    Create files for clusters containing 1, 2, 3, and 4 genes.
    """
    for gene_count in range(1, 5):
        filename = f"{gene_count}_gene_cluster.txt"
        clusters = [str(cluster) for cluster, count in cluster_counts.items() if count == gene_count]
        if clusters:
            with open(filename, 'w') as f:
                f.write(','.join(clusters))
            print(f"Created {filename} with {len(clusters)} clusters")
        else:
            print(f"No clusters found with exactly {gene_count} gene(s)")

def print_results(cluster_counts, total_genes):
    """
    Print the analysis results including distribution of genes per cluster.
    """
    if not cluster_counts:
        print("\nNo matching genes found in any cluster.")
        return
        
    # Basic cluster results
    print("\nCluster Results:")
    print("-" * 50)
    print("Cluster | Matching genes | Percentage")
    print("-" * 50)
    
    total_matches = sum(cluster_counts.values())
    
    for cluster, count in sorted(cluster_counts.items(), key=lambda item: item[1], reverse=False):
        count = cluster_counts[cluster]
        percentage = (count / total_genes) * 100
        print(f"{cluster:7d} | {count:13d} | {percentage:6.2f}%")
    
    print("-" * 50)
    print(f"Total genes checked: {total_genes}")
    print(f"Total matches found: {total_matches}")
    print(f"Number of clusters with matches: {len(cluster_counts)}")
    
    # Distribution analysis
    print("\nDistribution Analysis:")
    print("-" * 50)
    print("Genes per | Number of  | Portion of")
    print("cluster   | clusters   | total (%)")
    print("-" * 50)
    
    # Count how many clusters have X genes
    distribution = Counter(cluster_counts.values())
    total_clusters = len(cluster_counts)
    
    for gene_count in sorted(distribution.keys()):
        cluster_count = distribution[gene_count]
        portion = (cluster_count / total_clusters) * 100
        print(f"{gene_count:9d} | {cluster_count:10d} | {portion:9.2f}%")
    
    # Calculate some statistics
    max_genes = max(cluster_counts.values())
    min_genes = min(cluster_counts.values())
    avg_genes = total_matches / total_clusters
    
    print("-" * 50)
    print(f"Statistics:")
    print(f"Maximum genes in a cluster: {max_genes}")
    print(f"Minimum genes in a cluster: {min_genes}")
    print(f"Average genes per cluster: {avg_genes:.2f}")

    # Create cluster files
    create_cluster_files(cluster_counts)

def main():
    parser = argparse.ArgumentParser(
        description="Analyze gene distribution across clusters"
    )
    
    parser.add_argument(
        "cluster_file",
        help="Path to file containing gene IDs and cluster numbers"
    )
    
    parser.add_argument(
        "gene_list",
        help="Path to file containing list of genes to check"
    )
    
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug information"
    )
    
    args = parser.parse_args()
    
    try:
        with open(args.gene_list, 'r') as f:
            total_genes = sum(1 for line in f if line.strip())
            
        results = analyze_gene_clusters(args.cluster_file, args.gene_list, args.debug)
        print_results(results, total_genes)
        
    except Exception as e:
        print(f"An unexpected error occurred: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()