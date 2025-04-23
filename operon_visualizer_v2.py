#!/usr/bin/env python3

import os
import sys
import pandas as pd
import argparse
import colorsys

def get_predominant_strand(operon_genes, gff_data):
    strand_counts = {1: 0, -1: 0}
    for _, gene in operon_genes.iterrows():
        gff_gene = gff_data[gff_data['locus_tag'] == gene['locus_tag_x']]
        if not gff_gene.empty:
            strand_counts[gff_gene.iloc[0]['strand']] += 1
    return 1 if strand_counts[1] >= strand_counts[-1] else -1

def parse_gff(file_path):
    print(f"Parsing GFF file: {file_path}")
    genes = []
    with open(file_path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) == 9 and parts[2] == 'CDS':
                try:
                    attributes = {}
                    for item in parts[8].split(';'):
                        if '=' in item:
                            key, value = item.split('=', 1)
                            attributes[key] = value
                    genes.append({
                        'seqid': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': 1 if parts[6] == '+' else -1,
                        'locus_tag': attributes.get('locus_tag', ''),
                        'product': attributes.get('product', '')
                    })
                except Exception as e:
                    print(f"Warning: Skipping malformed line in {file_path}: {str(e)}")
                    continue
    print(f"Parsed {len(genes)} genes from GFF file")
    return pd.DataFrame(genes)

def parse_operon_data(file_path):
    print(f"Parsing operon data file: {file_path}")
    try:
        data = pd.read_csv(file_path, sep='\t', encoding='utf-8')
    except UnicodeDecodeError:
        data = pd.read_csv(file_path, sep='\t', encoding='latin-1')
    print(f"Parsed {len(data)} rows from operon data file")
    return data

def generate_colors(n):
    HSV_tuples = [(x*1.0/n, 0.8, 0.9) for x in range(n)]  # Increased saturation to 0.8 and value to 0.9
    RGB_tuples = [colorsys.hsv_to_rgb(*x) for x in HSV_tuples]
    return ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) for r, g, b in RGB_tuples]

def generate_legend():
    legend = """
    <h3>Color Legend</h3>
    <div style='display: flex; flex-wrap: wrap; margin-bottom: 20px;'>
        <div style='margin: 5px;'><span style='background-color: #FAF0D7; padding: 5px;'>No annotation</span></div>
        <div style='margin: 5px;'><span style='background-color: #BFA784; padding: 5px;'>Core gene</span></div>
        <div style='margin: 5px;'><span style='background-color: #DD5341; padding: 5px;'>Highlighted gene</span></div>
        <div style='margin: 5px;'><span style='background-color: #0C8D90; padding: 5px;'>AFP11</span></div>
        <div style='margin: 5px;'><span style='background-color: #387654; padding: 5px;'>AFP12</span></div>
    </div>
    """
    return legend

def generate_html(genome_id, operons, gff_data, color_map, highlight_locus_tags, strand=None):
    print(f"Generating HTML for genome: {genome_id}")
    html_content = f"""
    <!-- <h2>Genome: {genome_id}</h2> -->
    <style>
    .body {{
        margin: 0;
        padding: 0;
        line-height: 1.2;
    }}
    .operon-container {{
        margin-bottom: 2px;  /* Increased margin between operons */
        position: relative;  /* Contains absolutely positioned children */
        overflow: hidden;  /* Ensures container expands to fit content */
    }}
    .operon {{
        position: relative;
        height: 40px;  /* Fixed height for operon */
    }}
    .gene-container {{
        position: absolute;
        top: 50%;
        transform: translateY(-50%);
    }}
    .gene-background, .gene-foreground {{
        height: 30px;
        position: absolute;
        left: 0;
    }}
    /* ... rest of your existing CSS ... */
    </style>
    """

    SCALE = 3 / 100  # 10 pixels per 100 nucleotides
    OUTLINE_WIDTH = 2  # Width of the outline in pixels

    cluster_count = 0
    for cluster_id, cluster_genes in operons.groupby('clust_ID_x'):
        try:
            # Determine predominant strand for this operon
            predominant_strand = get_predominant_strand(cluster_genes, gff_data)
            
            # Filter based on strand if specified
            if strand is not None:
                if (strand == 'plus' and predominant_strand == -1) or (strand == 'minus' and predominant_strand == 1):
                    continue

            cluster_count += 1
            cluster_genes = cluster_genes.sort_values('start_coord_x')
            
            #core_or_highlighted_genes = cluster_genes[
            #    (cluster_genes['core_whole'].notna()) | 
            #    (cluster_genes['locus_tag_x'].isin(highlight_locus_tags))
            #]
            
            #if core_or_highlighted_genes.empty:
            #    continue
            
            #start_index = cluster_genes.index.get_loc(core_or_highlighted_genes.index[0])
            #end_index = cluster_genes.index.get_loc(core_or_highlighted_genes.index[-1])
            #cluster_genes = cluster_genes.iloc[start_index:end_index+1]
            
            total_width = sum((gene['end_coord_x'] - gene['start_coord_x']) * SCALE for _, gene in cluster_genes.iterrows())

            html_content += f"<div class='operon-container'>"  # New container div
            html_content += f"<h3 style='margin-bottom: 5px;'>Cluster ID: {cluster_id}</h3>"
            html_content += f"<div class='operon' style='width: {total_width}px;'>"

            operon_start = cluster_genes['start_coord_x'].min()

            for _, gene in cluster_genes.iterrows():
                gene_length = gene['end_coord_x'] - gene['start_coord_x']
                gene_width = gene_length * SCALE
                gene_position = (gene['start_coord_x'] - operon_start) * SCALE
                
                gff_gene = gff_data[gff_data['locus_tag'] == gene['locus_tag_x']]
                
                # Updated color selection logic
                if gene['locus_tag_x'] in highlight_locus_tags:
                    color = '#DD5341'  # Highlighted gene color
                elif pd.notna(gene['Short name']) and 'Tail_P2_I' in str(gene['Short name']):
                    color = '#FFA500'  # Orange for Tail_P2_I genes
                elif pd.isna(gene['core_whole']):
                    color = '#FAF0D7'  # No annotation color
                elif gene['core_whole'] == 'afp11':
                    color = '#0C8D90'  # AFP11 color
                elif gene['core_whole'] == 'afp12':
                    color = '#387654'  # AFP12 color
                else:
                    color = '#BFA784'  # Core gene color
                
                if not gff_gene.empty:
                    gff_gene = gff_gene.iloc[0]
                    strand_class = 'forward' if gff_gene['strand'] == 1 else 'reverse'
                    html_content += f"""
                        <div class='gene-container' style='width: {gene_width + OUTLINE_WIDTH}px; left: {gene_position}px;'>
                            <div class='gene-background {strand_class}-background' 
                                 style='width: {gene_width + OUTLINE_WIDTH}px;'>
                            </div>
                            <div class='gene-foreground {strand_class}-foreground' 
                                 style='width: {gene_width}px; background-color: {color};' 
                                 title='{gff_gene['product']}'>
                            </div>
                        </div>
                    """
                else:
                    html_content += f"""
                        <div class='gene-container' style='width: {gene_width + OUTLINE_WIDTH}px; left: {gene_position}px;'>
                            <div class='gene-background {strand_class}-background' 
                                 style='width: {gene_width + OUTLINE_WIDTH}px;'>
                            </div>
                            <div class='gene-foreground {strand_class}-foreground' 
                                 style='width: {gene_width}px; background-color: {color};'>
                            </div>
                        </div>
                    """

            html_content += "</div></div>"  # Close both operon and operon-container divs

        except Exception as e:
            print(f"Warning: Skipping operon {cluster_id} due to error: {str(e)}")
            continue

    return html_content

def main(base_dir, operon_file, output_file, clust_ids=None, locus_tags_file=None, strand=None):


    print("Starting operon visualization process")
    
    operon_data = parse_operon_data(operon_file)
    
    # Filter operons based on specified cluster IDs
    if clust_ids:
        cluster_list = [int(cid.strip()) for cid in clust_ids.split(',')]
        operon_data = operon_data[operon_data['clust_ID_x'].isin(cluster_list)]
        print(f"Filtered operon data to include {len(cluster_list)} specified clusters")
    
    # Read locus tags file if provided
    highlight_locus_tags = set()
    if locus_tags_file:
        with open(locus_tags_file, 'r') as f:
            highlight_locus_tags = set(line.strip() for line in f)
        print(f"Loaded {len(highlight_locus_tags)} locus tags to highlight")

    # Filter operons based on specified cluster IDs
    if clust_ids:
        cluster_list = [int(cid.strip()) for cid in clust_ids.split(',')]
        operon_data = operon_data[operon_data['clust_ID_x'].isin(cluster_list)]
        print(f"Filtered operon data to include {len(cluster_list)} specified clusters")
    
    html_content = """
    <html>
    <head>
        <style>
        .operon {
            display: flex;
            align-items: center;
            margin-bottom: 20px;
        }
        .gene-container {
            position: absolute;
        }
        .gene-background, .gene-foreground {
            height: 30px;
            position: absolute;
            top: 0;
            left: 0;
        }
        .gene-background {
            background-color: #000; /* Color of the outline */
            z-index: 1;
        }
        .gene-foreground {
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            z-index: 2;
        }
        .forward-background {
            clip-path: polygon(0% 0%, 90% 0%, 100% 50%, 90% 100%, 0% 100%);
        }
        .forward-foreground {
            clip-path: polygon(1px 1px, calc(90% - 1px) 1px, calc(100% - 1px) 50%, calc(90% - 1px) calc(100% - 1px), 1px calc(100% - 1px));
        }
        .reverse-background {
            clip-path: polygon(10% 0%, 100% 0%, 100% 100%, 10% 100%, 0% 50%);
        }
        .reverse-foreground {
            clip-path: polygon(calc(10% + 1px) 1px, calc(100% - 1px) 1px, calc(100% - 1px) calc(100% - 1px), calc(10% + 1px) calc(100% - 1px), 1px 50%);
        }
        </style>
    </head>
    <body>
        <h1>Operon Visualizations</h1>
    """
    
    html_content += "<h2>Color Legend</h2>"
    html_content += "<div style='display: flex; flex-wrap: wrap;'>"
    html_content += "<div style='margin: 5px;'><span style='background-color: #FAF0D7; padding: 5px;'>No annotation</span></div>"
    html_content += "<div style='margin: 5px;'><span style='background-color: #BFA784; padding: 5px;'>Core gene</span></div>"
    html_content += "<div style='margin: 5px;'><span style='background-color: #811638; padding: 5px;'>Highlighted gene</span></div>"
    html_content += "<div style='margin: 5px;'><span style='background-color: #0C8D90; padding: 5px;'>AFP11</span></div>/span></div>"
    html_content += "<div style='margin: 5px;'><span style='background-color: #0B7978; padding: 5px;'>AFP12</span></div>/span></div>"
    html_content += "<div style='margin: 5px;'><span style='background-color: #FFA500; padding: 5px;'>Tail_P2_I gene</span></div>"


    html_content += "</div>"

    genome_count = 0
    for genome_id in operon_data['genome_ID_x'].unique():
        genome_count += 1
        print(f"Processing genome {genome_count}: {genome_id}")
        genome_dir = os.path.join(base_dir, str(genome_id))
        if os.path.isdir(genome_dir):
            gff_file = os.path.join(genome_dir, f"{genome_id}.gff")
            if os.path.exists(gff_file):
                gff_data = parse_gff(gff_file)
                genome_operons = operon_data[operon_data['genome_ID_x'] == genome_id]
                html_content += generate_html(genome_id, genome_operons, gff_data, None, highlight_locus_tags, strand)
            else:
                print(f"Warning: GFF file not found for genome {genome_id}")
        else:
            print(f"Warning: Directory not found for genome {genome_id}")

    html_content += """
    </body>
    </html>
    """

    with open(output_file, 'w') as f:
        f.write(html_content)

    print(f"Processed {genome_count} genomes")
    print(f"HTML file saved as {output_file}")
    print("Operon visualization process completed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate operon visualizations from GFF files and operon data.")
    parser.add_argument("input_dir", help="Directory containing genome folders with GFF files")
    parser.add_argument("operon_file", help="Tab-separated file containing operon data")
    parser.add_argument("output_file", help="Output HTML file for operon visualizations (will add .html if not present)")
    parser.add_argument("--clust", help="Comma-separated list of cluster IDs to include", default=None)
    parser.add_argument("--locus_tags", help="File containing locus tags to highlight", default=None)
    parser.add_argument("--strand", choices=['plus', 'minus'], help="Filter operons by predominant strand orientation")
    args = parser.parse_args()

    # Ensure the output file has a .html extension
    output_file = args.output_file if args.output_file.endswith('.html') else args.output_file + '.html'

    main(args.input_dir, args.operon_file, output_file, args.clust, args.locus_tags, args.strand)