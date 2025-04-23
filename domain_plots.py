import os
import random
import colorsys
import argparse
from collections import defaultdict
from Bio import SeqIO
from typing import List, Dict, Tuple
import glob

def parse_data_files(file_patterns: List[str], score_threshold: float) -> List[Tuple]:
    """
    Parse multiple CATH hits files using pattern matching.
    
    Args:
        file_patterns: List of file patterns (can include wildcards)
        score_threshold: Minimum score to include a domain
    
    Returns:
        List of parsed records
    """
    records = []
    processed_files = set()
    
    # Process each pattern
    for pattern in file_patterns:
        # Expand pattern to get matching files
        matching_files = glob.glob(pattern)
        if not matching_files:
            print(f"Warning: No files found matching pattern '{pattern}'")
            continue
            
        for file_path in matching_files:
            if file_path in processed_files:
                continue
                
            if not os.path.exists(file_path):
                print(f"Warning: File {file_path} does not exist, skipping...")
                continue
            
            print(f"Processing file: {file_path}")
            try:
                with open(file_path, 'r') as f:
                    lines = f.readlines()

                for line in lines:
                    line = line.strip()
                    fields = line.split()
                    if len(fields) == 7:
                        gene_id = fields[0]
                        domain_id = fields[1]
                        score = float(fields[2])
                        if score < score_threshold:
                            continue
                        hit_coord1 = fields[3]
                        hit_coord2 = fields[4]
                        evalue1 = float(fields[5])
                        evalue2 = float(fields[6])
                        records.append((gene_id, domain_id, score, hit_coord1, hit_coord2, evalue1, evalue2))
                
                processed_files.add(file_path)
                
            except Exception as e:
                print(f"Error processing file {file_path}: {str(e)}")
                continue
    
    print(f"Successfully processed {len(processed_files)} files")
    print(f"Total records found: {len(records)}")
    return records

def read_fasta(fasta_file):
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"The file {fasta_file} does not exist.")
    
    protein_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id.split()[0]
        protein_lengths[gene_id] = len(record.seq)
    return protein_lengths

def resolve_overlaps(domain_records, custom_overlap_threshold=0.6):
    resolved_records = []
    custom_domains = defaultdict(list)
    pfam_domains = defaultdict(list)

    for record in domain_records:
        gene_id, domain_id, _, hit_coord1, _, _, _ = record
        if domain_id.startswith("domain"):
            custom_domains[gene_id].append(record)
        else:
            pfam_domains[gene_id].append(record)

    for gene_id in custom_domains:
        for custom_record in custom_domains[gene_id]:
            _, custom_domain_id, _, custom_hit_coord1, _, _, _ = custom_record
            overlap_resolved = False
            for pfam_record in pfam_domains[gene_id]:
                _, pfam_domain_id, _, pfam_hit_coord1, _, _, _ = pfam_record
                for custom_coord in custom_hit_coord1.split(","):
                    custom_start, custom_end = map(int, custom_coord.split("-"))
                    custom_length = custom_end - custom_start + 1
                    for pfam_coord in pfam_hit_coord1.split(","):
                        pfam_start, pfam_end = map(int, pfam_coord.split("-"))

                        overlap_start = max(custom_start, pfam_start)
                        overlap_end = min(custom_end, pfam_end)
                        overlap_length = max(0, overlap_end - overlap_start + 1)

                        if overlap_length / custom_length >= custom_overlap_threshold:
                            resolved_records.append((gene_id, f"{pfam_domain_id}({custom_domain_id})", *custom_record[2:]))
                            overlap_resolved = True
                            break
                    if overlap_resolved:
                        break
                if overlap_resolved:
                    break
            if not overlap_resolved:
                resolved_records.append(custom_record)

    for gene_id in pfam_domains:
        for pfam_record in pfam_domains[gene_id]:
            resolved_records.append(pfam_record)

    return resolved_records

def get_palette_color():
    """Generate a color from the combined palette."""
    color_palette = [
        "#fee327", "#fdca54", "#f6a570", 
        "#f1969b", "#f08ab1", "#c78dbd", 
        "#927db6", "#5da0d7", "#00b3e1", 
        "#50bcbf", "#65bda5", "#87bf54",
        "#881177", "#aa3355", "#cc6666",
        "#ee9944", "#eedd00", "#99dd55",
        "#44dd88", "#22ccbb", "#00bbcc",
        "#0099cc", "#3366bb", "#663399"
    ]
    return random.choice(color_palette)

def generate_html(records, protein_lengths, fasta_file, output_file, presentation_scale=0.8):
    def get_text_color(hex_color):
        hex_color = hex_color.lstrip('#')
        r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
        brightness = (r * 299 + g * 587 + b * 114) / 1000
        return '#FFFFFF' if brightness < 128 else '#000000'

    def get_outline_color(text_color):
        r, g, b = int(text_color[1:3], 16), int(text_color[3:5], 16), int(text_color[5:7], 16)
        outline_r = 255 - r
        outline_g = 255 - g
        outline_b = 255 - b
        return f"rgb({outline_r}, {outline_g}, {outline_b})"

    colors = {}
    domain_names = set(record[1] for record in records)
    for domain in domain_names:
        colors[domain] = get_palette_color()

    organism_data = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                parts = line.strip().split()
                gene_id = parts[0][1:]
                full_description = " ".join(parts[1:])
                organism_data[gene_id] = full_description

    gene_groups = {}
    for record in records:
        gene_id, domain_id, _, _, _, _, _ = record
        if gene_id not in gene_groups:
            gene_groups[gene_id] = domain_id

    sorted_gene_groups = sorted(gene_groups.items(), key=lambda x: int(records[next(i for i, rec in enumerate(records) if rec[0] == x[0])][3].split("-")[0]))

    displayed_genes = set()

    with open(output_file, 'w') as f:
        # HTML header and style definitions
        f.write("""<!DOCTYPE html>
<html lang='en'>
<head>
<meta charset='UTF-8'>
<title>Protein Domain Architecture</title>
<link rel='stylesheet' href='https://cdnjs.cloudflare.com/ajax/libs/tippy.js/6.3.1/tippy.css' />
<style>
.protein { margin: 20px 0; font-family: Arial, sans-serif; font-size: 14px; white-space: nowrap; position: relative; }
.domain { border: 2px solid #8BA0A4; display: inline-flex; justify-content: center; align-items: center; height: 20px; margin: 0; position: relative; top: 5px; border-radius: 5px; box-shadow: 0 0 1px rgba(0,0,0,0.1); }
.domain-label { font-size: 12px; font-weight: bold;}
.tooltip { font-size: 14px; color: #fff; background-color: #333; padding: 5px; border-radius: 3px; }
.gap { display: inline-block; height: 2px; background-color: #8BA0A4; position: relative; top: 4px; margin: 0; transform: translateY(-50%); z-index: -1; }
.gene-id { font-family: Arial, sans-serif; font-size: 16px; font-weight: bold; }
.domain-text { font-weight: bold; position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); z-index: 1; background-color: inherit; padding: 1px 2px; border-radius: 4px;}
</style>
</head>
<body>
<script src='https://cdnjs.cloudflare.com/ajax/libs/tippy.js/6.3.1/tippy.umd.min.js'></script>
""")

        for gene_id, domain_id in sorted_gene_groups:
            if gene_id in displayed_genes:
                continue

            gene_ids = [gene for gene, dom in gene_groups.items() if dom == domain_id]
            f.write(f"<h2>Genes with First Domain: {domain_id}</h2>\n")

            subgroups = defaultdict(list)
            for gene_id in gene_ids:
                subgroup_key = tuple(rec[1] for rec in records if rec[0] == gene_id)
                subgroups[subgroup_key].append(gene_id)

            for subgroup_key, subgroup_gene_ids in subgroups.items():
                f.write(f"<details open><summary>Domain Architecture: {' - '.join(subgroup_key)} ({len(subgroup_gene_ids)} genes)</summary>\n")
                f.write("<div id='protein-container'>\n")
                
                for gene_id in subgroup_gene_ids:
                    displayed_genes.add(gene_id)
                    full_description = organism_data.get(gene_id, "Unknown")
                    f.write(f"<div><span class='gene-id'>{gene_id}</span> {full_description}</div>\n")
                    f.write("<div class='protein'>\n")

                    domain_records = [record for record in records if record[0] == gene_id]
                    last_end = 0
                    protein_length = protein_lengths.get(gene_id, last_end)
                    
                    for _, domain_id, _, hit_coord1, _, _, _ in domain_records:
                        for coord in hit_coord1.split(","):
                            start, end = map(int, coord.split("-"))
                            domain_length = (end - start + 1) * presentation_scale
                            if start > last_end:
                                gap_length = (start - last_end) * presentation_scale
                                f.write(f"<span class='gap' style='width: {gap_length}px;'></span>\n")
                            text_color = get_text_color(colors[domain_id])
                            outline_color = get_outline_color(text_color)
                            f.write(f"<span class='domain' style='width: {domain_length}px; background-color: {colors[domain_id]}; position: relative;'>")
                            f.write(f"<span class='domain-text' style='color: {text_color};'>{domain_id}</span>")
                            f.write(f"<span class='domain-outline' style='color: {outline_color}; text-shadow: 1px 1px 0 {outline_color}, -1px -1px 0 {outline_color}, 1px -1px 0 {outline_color}, -1px 1px 0 {outline_color};'>{domain_id}</span></span>\n")
                            last_end = end

                    if last_end < protein_length:
                        gap_length = (protein_length - last_end) * presentation_scale
                        f.write(f"<span class='gap' style='width: {gap_length}px;'></span>\n")

                    f.write(f"<span style='display: inline-block; vertical-align: middle; margin-left: 5px;'>({protein_length}aa)</span>")
                    f.write("</div>\n")
                f.write("</div>\n")
                f.write("</details>\n")

        f.write("<script>tippy('.domain');</script>\n")
        f.write("</body>\n</html>\n")

def main():
    parser = argparse.ArgumentParser(description="Generate HTML visualization of protein domain architecture from multiple CATH hits files.")
    parser.add_argument("--cath_hits", nargs='+', required=True, 
                      help="One or more CATH hits files or patterns (e.g., 'hits*.txt')")
    parser.add_argument("--faa", required=True, 
                      help="Path to the protein FASTA file")
    parser.add_argument("--output", required=True, 
                      help="Path for the output HTML file")
    parser.add_argument("--score_threshold", type=float, default=0,
                      help="Minimum score threshold for domain inclusion (default: 0)")
    parser.add_argument("--scale", type=float, default=0.8,
                      help="Presentation scale for domain visualization (default: 0.8)")
    
    args = parser.parse_args()

    # Parse all domain files
    records = parse_data_files(args.cath_hits, args.score_threshold)
    
    if not records:
        print("No valid records found in input files. Exiting...")
        return

    # Resolve domain overlaps
    resolved_records = resolve_overlaps(records)

    # Read protein lengths from FASTA
    protein_lengths = read_fasta(args.faa)

    # Generate HTML visualization
    generate_html(resolved_records, protein_lengths, args.faa, args.output, args.scale)
    print(f"HTML visualization has been saved to: {args.output}")

if __name__ == "__main__":
    main()