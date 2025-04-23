#!/usr/bin/env python3
import argparse
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from collections import defaultdict

def parse_protein_lengths(file_path):
    """Parse a tab-delimited file with protein IDs and lengths."""
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
    """Parse the cluster data file."""
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
                            'members': []
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

def extract_protein_base_id(member_id):
    """Extract the base protein ID from a member ID string."""
    match = re.match(r'(.+)_\d+_\d+(?:\s*$|\[)', member_id)
    if match:
        return match.group(1)
    return None

def calculate_domain_position(member_id):
    """Extract coordinates and calculate domain position."""
    match = re.search(r'(.+)_(\d+)_(\d+)(?:\s*$|\[)', member_id)
    if match:
        protein_id = match.group(1)
        coord1 = int(match.group(2))
        coord2 = int(match.group(3))
        
        # Ensure coord2 is greater than coord1
        if coord2 < coord1:
            coord1, coord2 = coord2, coord1
        
        domain_size = coord2 - coord1
        avg_position = (coord1 + coord2) / 2.0
        
        return protein_id, avg_position, domain_size
    
    return None, None, None

def find_protein_id_variations(protein_id, protein_lengths):
    """Try different variations of the protein ID to find a match."""
    # 1. Try exact match
    if protein_id in protein_lengths:
        return protein_id, protein_lengths[protein_id]
    
    # 2. Try removing last underscore part if present
    if '_' in protein_id:
        parts = protein_id.split('_')
        potential_id = '_'.join(parts[:-1])
        if potential_id in protein_lengths:
            return potential_id, protein_lengths[potential_id]
    
    # 3. Try just the first part (before any underscore)
    first_part = protein_id.split('_')[0]
    if first_part in protein_lengths:
        return first_part, protein_lengths[first_part]
    
    return None, None

def build_network_data(clusters, protein_lengths, max_coverage=0.6, domain_size_limit=600):
    """Build network data from cluster information."""
    nodes = {}
    protein_to_clusters = defaultdict(list)
    filtered_clusters = []
    
    # Process each cluster
    for idx, cluster in enumerate(clusters):
        domain_sizes = []
        avg_coords = []
        protein_sizes = []
        coords_with_sizes = []
        protein_base_ids = set()  # Use a set to avoid duplicates
        
        for member in cluster['members']:
            protein_id, avg_position, domain_size = calculate_domain_position(member)
            
            if protein_id is None or avg_position is None or domain_size is None:
                continue
            
            # Extract the base protein ID (for connecting with other clusters)
            protein_base_id = extract_protein_base_id(member)
            if protein_base_id:
                protein_base_ids.add(protein_base_id)
            
            # Try to find the protein length
            matched_id, protein_length = find_protein_id_variations(protein_id, protein_lengths)
            
            # Skip if domain is too large relative to protein length, or unreasonably large
            if matched_id is not None:
                if domain_size / protein_length > max_coverage:
                    continue
                protein_sizes.append(protein_length)
                coords_with_sizes.append(avg_position)
            elif domain_size > domain_size_limit:
                continue
            
            domain_sizes.append(domain_size)
            avg_coords.append(avg_position)
        
        # Skip clusters with no valid members
        if not domain_sizes or not avg_coords:
            continue
        
        # Calculate average values for this cluster
        avg_domain_size = np.mean(domain_sizes)
        overall_avg_coord = np.mean(avg_coords)
        
        # Calculate relative position if protein sizes are available
        if protein_sizes and coords_with_sizes:
            avg_protein_size = np.mean(protein_sizes)
            avg_coord_with_size = np.mean(coords_with_sizes)
            relative_position = min(1.0, max(0.0, avg_coord_with_size / avg_protein_size))
        else:
            avg_protein_size = 0
            relative_position = 0.5  # Default to middle
        
        nodes[idx] = {
            'size': cluster['size'],
            'avg_domain_size': avg_domain_size,
            'avg_coordinate': overall_avg_coord,
            'avg_protein_size': avg_protein_size,
            'color': relative_position,
            'label': f"{cluster['representative']}",
            'protein_ids': protein_base_ids
        }
        
        filtered_clusters.append(cluster)
        
        # Track which clusters each protein belongs to
        for protein_id in protein_base_ids:
            protein_to_clusters[protein_id].append(idx)
    
    # Build edges between clusters that share proteins
    edges = {}
    for protein_id, cluster_indices in protein_to_clusters.items():
        # This protein might belong to multiple clusters (has multiple domains)
        if len(cluster_indices) > 1:
            for i in range(len(cluster_indices)):
                for j in range(i+1, len(cluster_indices)):
                    cluster_i = cluster_indices[i]
                    cluster_j = cluster_indices[j]
                    
                    if cluster_i in nodes and cluster_j in nodes:  # Make sure both nodes exist
                        edge_key = tuple(sorted([cluster_i, cluster_j]))
                        if edge_key not in edges:
                            edges[edge_key] = {'weight': 0, 'shared_proteins': set()}
                        edges[edge_key]['weight'] += 1
                        edges[edge_key]['shared_proteins'].add(protein_id)
    
    return nodes, edges, protein_to_clusters, filtered_clusters

def create_network_plot(nodes, edges, output_file=None, min_edge_weight=1, label_edges=False, node_scale=1.0):
    """Create a network plot using NetworkX and Matplotlib."""
    # Create a new graph
    G = nx.Graph()
    
    # Add nodes with attributes
    for node_id, attrs in nodes.items():
        G.add_node(node_id, 
                   size=attrs['size'],
                   color=attrs['color'],
                   label=attrs['label'])
    
    # Add edges with weights (filter by minimum weight)
    for (i, j), attrs in edges.items():
        if len(attrs['shared_proteins']) >= min_edge_weight:
            G.add_edge(i, j, weight=len(attrs['shared_proteins']), shared_proteins=len(attrs['shared_proteins']))
    
    # Remove nodes with no connections
    isolated_nodes = list(nx.isolates(G))
    G.remove_nodes_from(isolated_nodes)
    if isolated_nodes:
        print(f"Removed {len(isolated_nodes)} isolated nodes", file=sys.stderr)
    
    if len(G.nodes()) == 0:
        print("No nodes left after filtering. Try lowering min_edge_weight.", file=sys.stderr)
        return
    
    # Set up the plot
    plt.figure(figsize=(10, 7))
    
    pos = nx.spring_layout(G, iterations=50, k=0.25)
    
    # Prepare node attributes for visualization
    node_sizes = [G.nodes[n]['size'] * 20 * node_scale for n in G.nodes()]
    node_colors = [G.nodes[n]['color'] for n in G.nodes()]
    
    # Prepare edge attributes for visualization
    edge_weights = [G.edges[e]['weight'] / 4 for e in G.edges()]
    
    # Draw the network
    nodes_drawn = nx.draw_networkx_nodes(G, pos, 
                                         node_size=node_sizes,
                                         node_color=node_colors, 
                                         cmap=plt.cm.rainbow,
                                         alpha=0.8)
    
    edges_drawn = nx.draw_networkx_edges(G, pos, 
                                         width=edge_weights,
                                         alpha=0.5,
                                         edge_color='gray')
    
    # Get the top 6-7 largest nodes
    node_sizes = [(n, G.nodes[n]['size']) for n in G.nodes()]
    node_sizes.sort(key=lambda x: x[1], reverse=True)  # Sort by size in descending order
    top_nodes = [n for n, _ in node_sizes[:7]]  # Take the top 7 nodes

    if top_nodes:
        nx.draw_networkx_labels(G, pos, 
                                labels={n: G.nodes[n]['label'] for n in top_nodes},
                                font_size=8,
                                font_weight='bold')
    
    # Optional: Draw edge labels (number of shared proteins)
    if label_edges:
        edge_labels = {(i, j): G.edges[(i, j)]['shared_proteins'] for i, j in G.edges()}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)
    
    # Add a color bar
    cbar = plt.colorbar(nodes_drawn, label="Relative Position on Protein (N to C-terminal)")
    
    # Add a size legend
    size_range = [500,2000,5000]  # Will display as 1, 20, 100 after division by 5
    legend_elements = []
    for size in size_range:
        legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                          markerfacecolor='gray', markersize=np.sqrt(size/5) * np.sqrt(node_scale),
                                          label=f'{size/5:.0f}'))
    
    # Add the legend
    plt.legend(handles=legend_elements, title="Cluster Size", loc="upper right")
    
    # Add title and adjust layout
    plt.title("Protein Domain Network\nNodes: Clusters, Edges: Shared Proteins")
    plt.axis('off')
    plt.tight_layout()
    
    # Save the figure
    if output_file:
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        print(f"Network plot saved to {output_file}")
    
    plt.show()

def main():
    parser = argparse.ArgumentParser(
        description="Generate a network plot from protein clustering data."
    )
    parser.add_argument("input_file", help="Path to the input data file")
    parser.add_argument("--output", help="Output file path for the network plot (e.g. network.png)", default="network.png")
    parser.add_argument("--lengths", help="Path to a TSV file with protein IDs and lengths", required=True)
    parser.add_argument("--max-coverage", type=float, default=0.6, 
                        help="Maximum allowed domain coverage (domain size / protein length). Default: 0.6")
    parser.add_argument("--domain-size-limit", type=int, default=600,
                        help="Maximum allowed domain size for proteins without length information. Default: 600")
    parser.add_argument("--min-edge-weight", type=int, default=1,
                        help="Minimum number of shared proteins to create an edge. Default: 1")
    parser.add_argument("--label-edges", action="store_true",
                        help="Label edges with the number of shared proteins")
    parser.add_argument("--node-scale", type=float, default=1.0,
                        help="Scale factor for node sizes (1.0 = 100%%). Default: 1.0")
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
    
    print(f"Found {len(clusters)} clusters. Building network data...", file=sys.stderr)
    nodes, edges, protein_to_clusters, filtered_clusters = build_network_data(
        clusters, 
        protein_lengths, 
        max_coverage=args.max_coverage,
        domain_size_limit=args.domain_size_limit
    )
    
    print(f"Network built with {len(nodes)} nodes and {len(edges)} edges.", file=sys.stderr)
    print(f"Found {len(protein_to_clusters)} unique proteins across all clusters.", file=sys.stderr)
    
    print(f"Creating network plot...", file=sys.stderr)
    create_network_plot(nodes, edges, 
                       output_file=args.output,
                       min_edge_weight=args.min_edge_weight,
                       label_edges=args.label_edges,
                       node_scale=args.node_scale)

if __name__ == "__main__":
    main()
