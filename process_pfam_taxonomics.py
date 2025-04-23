import argparse
import os
import csv
from tqdm import tqdm

def load_pfam_dict(pfam_file):
    pfam_dict = {}
    with open(pfam_file, 'r') as f:
        lines = f.readlines()
        for line in tqdm(lines, desc="Loading Pfam dictionary"):
            pfam_id, pfam_name = line.strip().split(maxsplit=1)
            pfam_dict[pfam_id] = pfam_name
    return pfam_dict

def load_taxa_dict(taxa_file):
    taxa_dict = {}
    with open(taxa_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in tqdm(reader, desc="Loading Taxa dictionary"):
            taxa_dict[row['TaxID']] = row
    return taxa_dict

def process_taxonomy_file(file_path, pfam_dict, taxa_dict, output_writer):
    pfam_id = os.path.basename(file_path).split('_')[0]
    pfam_name = pfam_dict.get(pfam_id, "NaN")
    
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            gene_id, taxa_id = line.strip().split('\t')
            taxa_info = taxa_dict.get(taxa_id, {})
            
            row = [
                pfam_id,
                pfam_name,
                gene_id,
                taxa_id,
                taxa_info.get('superkingdom', 'NaN'),
                taxa_info.get('kingdom', 'NaN'),  # Added kingdom
                taxa_info.get('phylum', 'NaN'),
                taxa_info.get('class', 'NaN'),
                taxa_info.get('order', 'NaN'),
                taxa_info.get('family', 'NaN'),
                taxa_info.get('genus', 'NaN'),
                taxa_info.get('species', 'NaN')
            ]
            output_writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Process taxonomy files and generate TSV output")
    parser.add_argument("--inputdir", required=True, help="Directory containing taxonomy files")
    parser.add_argument("--pfamdic", required=True, help="Pfam dictionary file")
    parser.add_argument("--taxadic", required=True, help="Taxonomic dictionary file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--output_header", required=False, help="Custom output header (comma-separated)")
    args = parser.parse_args()

    print("Loading dictionaries...")
    pfam_dict = load_pfam_dict(args.pfamdic)
    taxa_dict = load_taxa_dict(args.taxadic)

    print("Processing taxonomy files...")
    with open(args.output, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        
        if args.output_header:
            header = args.output_header.split(',')
        else:
            header = ['pfamID', 'pfam_name', 'gene_ID', 'gene_taxa_ID', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        
        writer.writerow(header)

        taxonomy_files = [f for f in os.listdir(args.inputdir) if f.endswith("_taxonomy.txt")]
        for filename in tqdm(taxonomy_files, desc="Processing files"):
            file_path = os.path.join(args.inputdir, filename)
            process_taxonomy_file(file_path, pfam_dict, taxa_dict, writer)

    print(f"Output written to {args.output}")

if __name__ == "__main__":
    main()
