#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import colorsys
import re
from collections import defaultdict
IMPORTANT_DOMAINS = [
    "Sheath", "Tail_P2_I", "afp11", "afp12", "Baseplate", "PAAR", "DUF4255", "DUF4157"
]
AA_TO_NT = 3
def parse_cluster_ranges(range_string):
    if not range_string:
        return None
    cluster_ids = set()
    parts = [part.strip() for part in range_string.split(',')]
    for part in parts:
        if '-' in part and not part.startswith('-'):
            try:
                start, end = part.split('-', 1)
                start_id = int(start.strip())
                end_id = int(end.strip())
                for cluster_id in range(start_id, end_id + 1):
                    cluster_ids.add(cluster_id)
            except ValueError:
                print(f"Warning: Invalid range format '{part}', skipping")
                continue
        else:
            try:
                cluster_id = int(part.strip())
                cluster_ids.add(cluster_id)
            except ValueError:
                print(f"Warning: Invalid cluster ID '{part}', skipping")
                continue
    result = sorted(list(cluster_ids))
    print(f"Parsed cluster range '{range_string}' -> {len(result)} cluster IDs: {result[:10]}{'...' if len(result) > 10 else ''}")
    return result
def parse_marker_gene_directory(marker_dir):
    marker_genes = {}
    if not marker_dir or not os.path.exists(marker_dir):
        return marker_genes
    print(f"Parsing marker gene directory: {marker_dir}")
    for filename in os.listdir(marker_dir):
        if filename.endswith('.txt'):
            marker_name = filename[:-4]
            filepath = os.path.join(marker_dir, filename)
            try:
                with open(filepath, 'r') as f:
                    gene_ids = set(line.strip() for line in f if line.strip())
                marker_genes[marker_name] = gene_ids
                print(f"  Loaded {len(gene_ids)} gene IDs for marker '{marker_name}'")
            except Exception as e:
                print(f"  Warning: Could not read {filename}: {str(e)}")
    print(f"Total markers loaded: {len(marker_genes)}")
    return marker_genes
def parse_jackhmmer_domtblout_multicoord(file_path):
    print(f"Parsing jackhmmer domtblout file with coordinates: {file_path}")
    gene_domains_multi = defaultdict(list)
    domain_set = set()
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                parts = line.strip().split()
                if len(parts) >= 21:
                    target_name = parts[0]
                    query_name = parts[3]
                    try:
                        # Standard HMMER3 domtblout format: columns 19, 20 are target ali_from, ali_to
                        t_ali_from = int(parts[19])
                        t_ali_to = int(parts[20])
                        gene_domains_multi[target_name].append({
                            'domain': query_name,
                            'start': t_ali_from,
                            'end': t_ali_to
                        })
                    except (ValueError, IndexError):
                        # If parsing coordinates fails, store domain without coordinates
                        gene_domains_multi[target_name].append({
                            'domain': query_name,
                            'start': None,
                            'end': None
                        })
                    domain_set.add(query_name)
                elif len(parts) >= 4:
                    # Fallback for files without coordinate columns
                    target_name = parts[0]
                    query_name = parts[3]
                    gene_domains_multi[target_name].append({
                        'domain': query_name,
                        'start': None,
                        'end': None
                    })
                    domain_set.add(query_name)
        print(f"Parsed {len(gene_domains_multi)} gene-domain mappings with {len(domain_set)} unique domains (with coordinates)")
        return gene_domains_multi, list(domain_set)
    except Exception as e:
        print(f"Error parsing jackhmmer file: {str(e)}")
        return defaultdict(list), []
def generate_domain_colors(domains):
    BUILTIN_COLORS = {
        'Sheath': '#8B4513',
        'Tail_P2_I': '#FFA500',
        'afp11': '#0C8D90',
        'afp12': '#387654'
    }
    colors = {}
    for domain in domains:
        if domain in BUILTIN_COLORS:
            colors[domain] = BUILTIN_COLORS[domain]
    remaining_domains = [d for d in domains if d not in BUILTIN_COLORS]
    if remaining_domains:
        n = len(remaining_domains)
        for i, domain in enumerate(sorted(remaining_domains)):
            hue = (i / n) * 0.7 + 0.15
            rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
            hex_color = '#%02x%02x%02x' % (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))
            colors[domain] = hex_color
    return colors
def parse_combined_gff(file_path, clust_ids=None):
    print(f"Parsing combined GFF file: {file_path}")
    try:
        data = pd.read_csv(file_path, sep='\t', encoding='utf-8')
    except UnicodeDecodeError:
        data = pd.read_csv(file_path, sep='\t', encoding='latin-1')
    if clust_ids:
        data = data[data['clusterID'].isin(clust_ids)]
        actual_clusters = sorted(data['clusterID'].unique())
        print(f"Filtered data to include clusters from specified ranges")
        print(f"Found {len(actual_clusters)} clusters in data: {actual_clusters[:10]}{'...' if len(actual_clusters) > 10 else ''}")
        if len(actual_clusters) == 0:
            print("WARNING: No clusters found matching the specified ranges!")
    genes = []
    for _, row in data.iterrows():
        try:
            attributes = {}
            for item in str(row['attributes']).split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    attributes[key] = value
            genes.append({
                'cluster_id': row['clusterID'],
                'genome_id': row['genomeNumber'],
                'gene_id': str(row['geneID']),
                'seqid': row['seqid'],
                'start': int(row['start']),
                'end': int(row['end']),
                'strand': 1 if row['strand'] == '+' else -1,
                'locus_tag': attributes.get('locus_tag', ''),
                'product': attributes.get('product', '')
            })
        except Exception as e:
            print(f"Warning: Skipping malformed line: {str(e)}")
            continue
    df = pd.DataFrame(genes)
    print(f"Parsed {len(df)} genes from {len(df['cluster_id'].unique())} clusters across {len(df['genome_id'].unique())} genomes")
    return df
def classify_clusters_by_markers(gff_data, marker_genes):
    cluster_classifications = defaultdict(lambda: defaultdict(set))
    if not marker_genes:
        for (genome_id, cluster_id), group in gff_data.groupby(['genome_id', 'cluster_id']):
            cluster_classifications['No markers'][genome_id].add(cluster_id)
        return dict(cluster_classifications)
    for (genome_id, cluster_id), cluster_data in gff_data.groupby(['genome_id', 'cluster_id']):
        cluster_gene_ids = set(cluster_data['gene_id'].tolist())
        present_markers = []
        for marker_name, marker_gene_set in marker_genes.items():
            if cluster_gene_ids & marker_gene_set:
                present_markers.append(marker_name)
        if present_markers:
            classification = ', '.join(sorted(present_markers))
        else:
            classification = 'No markers'
        cluster_classifications[classification][genome_id].add(cluster_id)
    return dict(cluster_classifications)
def filter_clusters_by_target(cluster_classifications, target_groups):
    filtered = {}
    for classification, genome_clusters in cluster_classifications.items():
        present_markers = [x.strip() for x in classification.split(',')]
        if any(marker in target_groups for marker in present_markers):
            filtered[classification] = genome_clusters
    return filtered
def get_operon_group_key(cluster_df, gene_domains_multi, marker_genes, important_domains, marker_names):
    domain_counts = []
    gene_ids = set(cluster_df['gene_id'])
    for dom in important_domains:
        count = 0
        for gid in gene_ids:
            for annot in gene_domains_multi.get(gid, []):
                if annot['domain'] == dom:
                    count += 1
                    break
        domain_counts.append(count)
    marker_counts = []
    for marker in marker_names:
        marker_counts.append(sum(gid in marker_genes.get(marker, set()) for gid in gene_ids))
    return tuple(domain_counts + marker_counts)
def calculate_domain_coverage(gene_domains_list, gene_length_aa):
    """
    Calculate total coverage for each domain across the gene.
    Returns dictionary with domain -> total_coverage_aa
    """
    domain_coverage = defaultdict(int)
    
    for annot in gene_domains_list:
        domain_name = annot['domain']
        if annot.get('start') is not None and annot.get('end') is not None:
            # Calculate coverage for this domain hit
            start_aa = max(1, annot['start'])
            end_aa = min(gene_length_aa, annot['end'])
            coverage = max(0, end_aa - start_aa + 1)
            domain_coverage[domain_name] += coverage
        else:
            # If no coordinates, assume full gene coverage for this domain
            domain_coverage[domain_name] += gene_length_aa
    
    return domain_coverage
def get_dominant_domain(gene_domains_list, gene_length_aa):
    """
    Returns the domain with the largest total coverage across the gene.
    """
    if not gene_domains_list:
        return None
    
    domain_coverage = calculate_domain_coverage(gene_domains_list, gene_length_aa)
    
    if not domain_coverage:
        return gene_domains_list[0]['domain']  # Fallback to first domain
    
    # Return domain with maximum coverage
    dominant_domain = max(domain_coverage.items(), key=lambda x: x[1])
    return dominant_domain[0]
def generate_html_for_cluster_multidomain(
        cluster_id, genome_id, cluster_genes, gene_domains_multi, domain_colors, marker_genes, marker_colors, domain_level=0
    ):
    SCALE = 3 / 100
    OUTLINE_WIDTH = 2
    DEFAULT_COLOR = '#BFA784'
    DEFAULT_MARKER_COLOR = '#DD5341'
    cluster_genes = cluster_genes.sort_values('start')
    total_width = sum((gene['end'] - gene['start']) * SCALE for _, gene in cluster_genes.iterrows())
    html_content = f"""
    <div class='operon-container'>
        <h4 style='margin: 2px 0;'>Cluster {cluster_id} (Genome: {genome_id})</h4>
        <div class='operon' style='width: {total_width}px; position:relative; height: 35px;'>
    """
    operon_start = cluster_genes['start'].min()
    for _, gene in cluster_genes.iterrows():
        gene_length_nt = gene['end'] - gene['start'] + 1
        gene_length_aa = gene_length_nt // 3  # Approximate gene length in amino acids
        gene_width_px = gene_length_nt * SCALE
        gene_position_px = (gene['start'] - operon_start) * SCALE
        gene_id = str(gene['gene_id'])
        
        # Check if gene is a marker
        is_marker = False
        marker_type = None
        for marker_name, marker_gene_set in marker_genes.items():
            if gene_id in marker_gene_set:
                is_marker = True
                marker_type = marker_name
                break
        
        gene_domains_list = gene_domains_multi.get(gene_id, [])
        
        # NEW LOGIC: Determine arrow color based on domain with largest coverage
        if is_marker:
            color = marker_colors.get(marker_type.lower(), DEFAULT_MARKER_COLOR)
            title_suffix = f" (Marker: {marker_type})"
        elif gene_domains_list:
            dominant_domain = get_dominant_domain(gene_domains_list, gene_length_aa)
            color = domain_colors.get(dominant_domain, DEFAULT_COLOR)
            
            # Calculate coverage info for title
            domain_coverage = calculate_domain_coverage(gene_domains_list, gene_length_aa)
            coverage_info = []
            for domain_name, coverage in sorted(domain_coverage.items(), key=lambda x: x[1], reverse=True):
                coverage_pct = (coverage / gene_length_aa) * 100 if gene_length_aa > 0 else 0
                coverage_info.append(f"{domain_name}({coverage_pct:.1f}%)")
            
            title_suffix = f" (Domains by coverage: {', '.join(coverage_info)})"
        else:
            color = DEFAULT_COLOR
            title_suffix = ""
        
        strand_class = 'forward' if gene['strand'] == 1 else 'reverse'
        
        # Generate domain/marker bars if domain_level is enabled
        bar_html = ''
        if domain_level and (is_marker or gene_domains_list):
            # For marker genes: show bar for first 150 amino acids
            if is_marker:
                aa_len = 150
                marker_bar_len_nt = min(aa_len * AA_TO_NT, gene_length_nt)
                bar_width_px = (marker_bar_len_nt / gene_length_nt) * gene_width_px
                bar_html += (
                    f"<div class='domain-bold-line' title='Marker: {marker_type} (first 150aa)' "
                    f"style='background-color: {color}; height: 4px; position: absolute; top: -8px; "
                    f"left: 0px; width: {bar_width_px}px; border-radius: 2px; z-index: 10;'></div>"
                )
            
            # For domain annotations: show bars for each domain region
            for annot in gene_domains_list:
                dom_col = domain_colors.get(annot['domain'], DEFAULT_COLOR)
                if annot.get('start') is not None and annot.get('end') is not None:
                    dom_start_aa = annot['start']
                    dom_end_aa = annot['end']
                    # Convert protein coordinates to nucleotide coordinates
                    dom_start_nt = max(0, (dom_start_aa - 1) * AA_TO_NT)
                    dom_end_nt = min(gene_length_nt, dom_end_aa * AA_TO_NT)
                    dom_nt_len = dom_end_nt - dom_start_nt
                    
                    if dom_nt_len > 0:
                        # If marker gene, avoid overlap with marker bar
                        if is_marker:
                            marker_bar_len_nt = min(150 * AA_TO_NT, gene_length_nt)
                            if dom_start_nt < marker_bar_len_nt:
                                if dom_end_nt > marker_bar_len_nt:
                                    # Partial overlap: show domain bar after marker bar
                                    dom_start_nt = marker_bar_len_nt
                                    dom_nt_len = dom_end_nt - dom_start_nt
                                else:
                                    # Complete overlap: skip this domain bar
                                    continue
                        
                        dom_px_left = (dom_start_nt / gene_length_nt) * gene_width_px
                        dom_px_width = (dom_nt_len / gene_length_nt) * gene_width_px
                        
                        if dom_px_width > 0:
                            bar_html += (
                                f"<div class='domain-bold-line' title='Domain: {annot['domain']} [{annot['start']}-{annot['end']} aa]' "
                                f"style='background-color: {dom_col}; height: 4px; position: absolute; top: -8px; "
                                f"left: {dom_px_left}px; width: {dom_px_width}px; border-radius: 2px; z-index: 10;'></div>"
                            )
                else:
                    # No coordinates available: show bar for entire gene if not a marker
                    if not is_marker:
                        bar_html += (
                            f"<div class='domain-bold-line' title='Domain: {annot['domain']} (full gene)' "
                            f"style='background-color: {dom_col}; height: 4px; position: absolute; top: -8px; "
                            f"left: 0px; width: {gene_width_px}px; border-radius: 2px; z-index: 10;'></div>"
                        )
        
        html_content += f"""
            <div class='gene-container clickable-gene' 
                 style='width: {gene_width_px + OUTLINE_WIDTH}px; left: {gene_position_px}px; position: absolute;'
                 data-gene-id='{gene_id}'
                 onclick='copyGeneId(this)'
                 title='Gene ID: {gene["gene_id"]} | {gene["locus_tag"]}: {gene["product"]}{title_suffix} | Click to copy Gene ID'>
                {bar_html}
                <div class='gene-background {strand_class}-background' 
                     style='width: {gene_width_px + OUTLINE_WIDTH}px;'>
                </div>
                <div class='gene-foreground {strand_class}-foreground' 
                     style='width: {gene_width_px}px; background-color: {color};'>
                </div>
            </div>
        """
    html_content += "</div></div>"
    return html_content
def generate_full_html(gff_data, gene_domains_multi, domain_colors, marker_genes, marker_colors, cluster_classifications, operon_order, group_labels, domain_level=0):
    html_content = """
    <html>
    <head>
        <style>
        body {
            margin: 10px;
            padding: 0;
            line-height: 1.1;
            font-family: Arial, sans-serif;
        }
        .section {
            margin-bottom: 15px;
        }
        .operon-container {
            margin-bottom: 8px;
            position: relative;
            overflow: visible;
        }
        .operon {
            position: relative;
            height: 35px;
            margin-top: 10px;
        }
        .gene-container {
            position: absolute;
        }
        .clickable-gene {
            cursor: pointer;
            transition: opacity 0.2s;
        }
        .clickable-gene:hover {
            opacity: 0.8;
        }
        .clickable-gene:active {
            transform: scale(0.95);
        }
        .gene-background, .gene-foreground {
            height: 25px;
            position: absolute;
            top: 0;
            left: 0;
        }
        .gene-background {
            background-color: #000;
            z-index: 1;
        }
        .gene-foreground {
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            z-index: 2;
            font-size: 9px;
        }
        .domain-bold-line {
            height: 4px;
            position: absolute;
            top: -8px;
            border-radius: 2px;
            z-index: 10;
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
        h1 { margin-bottom: 10px; }
        h2 { margin: 5px 0; color: #2c3e50; }
        h3 { margin: 8px 0 5px 0; color: #34495e; }
        h4 { margin: 2px 0; color: #7f8c8d; font-size: 14px; }
        .copy-notification {
            position: fixed;
            top: 20px;
            right: 20px;
            background-color: #4CAF50;
            color: white;
            padding: 10px 20px;
            border-radius: 5px;
            z-index: 1000;
            opacity: 0;
            transition: opacity 0.3s;
        }
        .copy-notification.show {
            opacity: 1;
        }
        </style>
    </head>
    <body>
        <h1>Cluster Operon Visualizations</h1>
        <p><em>Click on any gene arrow to copy its Gene ID to clipboard</em></p>
        <p><em>Gene arrow colors are determined by: 1) Marker genes (red), 2) Domain with largest coverage across gene, 3) Default (beige)</em></p>
    """
    
    if domain_level:
        html_content += "<p><em>Bold bars above genes indicate marker regions (first 150aa) or domain-specific regions with coordinates from jackhmmer results.</em></p>"
    
    html_content += "<div id='copyNotification' class='copy-notification'>Gene ID copied to clipboard!</div>"
    
    if marker_genes:
        html_content += "<h2>Marker Gene Legend</h2>"
        html_content += "<div style='display: flex; flex-wrap: wrap; margin-bottom: 10px;'>"
        for marker_name in sorted(marker_genes.keys()):
            color = marker_colors.get(marker_name.lower(), '#DD5341')
            html_content += f"<div style='margin: 3px;'><span style='background-color: {color}; padding: 3px; color: white; font-weight: bold;'>{marker_name}</span></div>"
        html_content += "</div>"
    
    if domain_colors:
        html_content += "<h2>Domain Legend</h2>"
        html_content += "<div style='display: flex; flex-wrap: wrap; margin-bottom: 10px;'>"
        builtin_domains = ['Sheath', 'Tail_P2_I', 'afp11', 'afp12']
        for domain in builtin_domains:
            if domain in domain_colors:
                color = domain_colors[domain]
                html_content += f"<div style='margin: 3px;'><span style='background-color: {color}; padding: 3px; color: white; font-weight: bold;'>{domain}</span></div>"
        for domain, color in sorted(domain_colors.items()):
            if domain not in builtin_domains:
                html_content += f"<div style='margin: 3px;'><span style='background-color: {color}; padding: 3px; color: white; font-weight: bold;'>{domain}</span></div>"
        html_content += f"<div style='margin: 3px;'><span style='background-color: #BFA784; padding: 3px; color: white;'>No domain</span></div>"
        html_content += "</div>"
    
    total_clusters = 0
    last_label = None
    for group_label, (group_key, genome_id, cluster_id, cluster_df) in zip(group_labels, operon_order):
        if group_label != last_label:
            html_content += f"<div class='section'><h2>Clusters with {group_label}</h2>"
            last_label = group_label
        html_content += generate_html_for_cluster_multidomain(
            cluster_id, genome_id, cluster_df, gene_domains_multi, domain_colors, marker_genes, marker_colors, domain_level=domain_level
        )
        total_clusters += 1
    
    html_content += """
        <script>
        async function copyGeneId(element) {
            const geneId = element.getAttribute('data-gene-id');
            try {
                await navigator.clipboard.writeText(geneId);
                showCopyNotification();
                console.log('Gene ID copied to clipboard:', geneId);
            } catch (err) {
                fallbackCopyTextToClipboard(geneId);
                console.error('Failed to copy using Clipboard API, using fallback:', err);
            }
        }
        function fallbackCopyTextToClipboard(text) {
            const textArea = document.createElement("textarea");
            textArea.value = text;
            textArea.style.top = "0";
            textArea.style.left = "0";
            textArea.style.position = "fixed";
            document.body.appendChild(textArea);
            textArea.focus();
            textArea.select();
            try {
                const successful = document.execCommand('copy');
                if (successful) {
                    showCopyNotification();
                }
            } catch (err) {
                console.error('Fallback: Oops, unable to copy', err);
            }
            document.body.removeChild(textArea);
        }
        function showCopyNotification() {
            const notification = document.getElementById('copyNotification');
            notification.classList.add('show');
            setTimeout(() => {
                notification.classList.remove('show');
            }, 2000);
        }
        </script>
    </body>
    </html>
    """
    print(f"Generated visualizations for {total_clusters} clusters (ordered by annotation/target copy numbers within each group)")
    return html_content
def parse_marker_color_option(marker_color_arg):
    """
    Parses the --marker_color argument string which has format like:
    ebap1:#000000,ebap2:#002121
    Returns a dict with lowercase marker names as keys and color strings as values.
    Case insensitive on keys.
    """
    color_map = {}
    if not marker_color_arg:
        return color_map
    parts = marker_color_arg.split(',')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        match = re.match(r"^([^:]+):(#?[0-9a-fA-F]{6})$", part, re.IGNORECASE)
        if match:
            marker_name = match.group(1).lower()
            color = match.group(2).upper()
            if not color.startswith('#'):
                color = '#' + color
            color_map[marker_name] = color
        else:
            print(f"Warning: Invalid marker_color format '{part}', expected name:#rrggbb")
    return color_map
def main(gff_file, output_file, clust_range=None, jackhmmer_file=None, marker_dir=None, target=None, domain_level=0, marker_color=None):
    print("Starting cluster operon visualization process")
    clust_ids = parse_cluster_ranges(clust_range)
    marker_genes = parse_marker_gene_directory(marker_dir)
    marker_colors = parse_marker_color_option(marker_color)
    gene_domains_multi = defaultdict(list)
    domain_colors = {}
    
    if jackhmmer_file:
        gene_domains_multi, domains = parse_jackhmmer_domtblout_multicoord(jackhmmer_file)
        domain_colors = generate_domain_colors(domains)
        print(f"Generated colors for {len(domain_colors)} domains")
    
    gff_data = parse_combined_gff(gff_file, clust_ids)
    if gff_data.empty:
        print("No data found matching the specified criteria")
        return
    
    if marker_genes:
        total_matches = 0
        for marker_name, marker_gene_set in marker_genes.items():
            matches = set(gff_data['gene_id'].tolist()) & marker_gene_set
            print(f"Found {len(matches)} matches for marker '{marker_name}'")
            total_matches += len(matches)
        print(f"Total marker gene matches: {total_matches}")
    
    cluster_classifications = classify_clusters_by_markers(gff_data, marker_genes)
    target_groups = None
    if target:
        target_groups = [x.strip() for x in target.split(',') if x.strip()]
        cluster_classifications = filter_clusters_by_target(cluster_classifications, target_groups)
        print(f"Filtered to target group(s): {', '.join(target_groups)}")
        if not cluster_classifications:
            print("No clusters found for specified target group(s).")
            return
    
    marker_names = sorted(marker_genes.keys())
    operon_order = []
    group_labels = []
    for classification in sorted(cluster_classifications.keys()):
        genome_clusters = cluster_classifications[classification]
        group_operons = []
        for genome_id in sorted(genome_clusters.keys()):
            cluster_ids = sorted(genome_clusters[genome_id])
            for cluster_id in cluster_ids:
                cluster_df = gff_data[(gff_data['genome_id'] == genome_id) & (gff_data['cluster_id'] == cluster_id)]
                group_key = get_operon_group_key(cluster_df, gene_domains_multi, marker_genes, IMPORTANT_DOMAINS, marker_names)
                group_operons.append((group_key, genome_id, cluster_id, cluster_df))
        group_operons.sort(key=lambda x: (x[0], x[1], x[2]))
        operon_order.extend(group_operons)
        group_labels.extend([classification] * len(group_operons))
    
    html_content = generate_full_html(
        gff_data, gene_domains_multi, domain_colors, marker_genes, marker_colors, cluster_classifications, operon_order, group_labels, domain_level=domain_level
    )
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    print(f"HTML file saved as {output_file}")
    print("Visualization process completed")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate operon visualizations with domain coverage-based gene coloring and multi-domain gene annotation ordering.")
    parser.add_argument("gff_file", help="GFF-style file containing operonic gene data")
    parser.add_argument("output_file", help="Output HTML file for operon visualizations")
    parser.add_argument("--clust", help="Cluster ranges (e.g., '1-10,20-50,100')", default=None)
    parser.add_argument("--jackhmmer", help="Jackhmmer domtblout file for domain-based coloring", default=None)
    parser.add_argument("--markers", help="Directory containing marker gene files (*.txt)", default=None)
    parser.add_argument("--target", help="Comma-separated list of marker gene groups to include (e.g. 'eBAP1,eBAP2')", default=None)
    parser.add_argument("--domain_level", type=int, choices=[0,1], default=0, help="Show bold marker/domain-colored bars above annotated gene regions: 0=off (default), 1=on")
    parser.add_argument("--marker_color", help="Comma-separated marker colors, e.g., ebap1:#000000,ebap2:#002121", default=None)
    args = parser.parse_args()
    output_file = args.output_file if args.output_file.endswith('.html') else args.output_file + '.html'
    main(args.gff_file, output_file, args.clust, args.jackhmmer, args.markers, args.target, args.domain_level, args.marker_color)
