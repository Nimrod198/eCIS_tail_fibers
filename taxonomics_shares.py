import argparse
import sys
from collections import defaultdict

def analyze_tsv(input_file, output_file):
    pfam_data = defaultdict(lambda: defaultdict(int))
    pfam_names = {}
    total_counts = defaultdict(int)

    with open(input_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            pfam_id, pfam_name, superkingdom, kingdom = fields[0], fields[1], fields[4], fields[5]
            
            # Store pfam names for later use
            pfam_names[pfam_id] = pfam_name
            
            pfam_data[pfam_id]['total'] += 1
            total_counts[pfam_id] += 1
            
            if superkingdom == 'Eukaryota':
                pfam_data[pfam_id][f'Eukaryota-{kingdom}'] += 1
            else:
                pfam_data[pfam_id][superkingdom] += 1

    results = []
    for pfam_id, data in pfam_data.items():
        taxa_info = []
        for taxa, count in sorted(data.items(), key=lambda x: x[1], reverse=True):
            if taxa != 'total':
                percentage = (count / data['total']) * 100
                taxa_info.append(f"{percentage:.1f}%,{taxa}")
        
        # Get the correct pfam name from the dictionary
        pfam_name = pfam_names.get(pfam_id, "Unknown")
        results.append(f"{pfam_id}\t{pfam_name}\t{{{';'.join(taxa_info)}}}")

    with open(output_file, 'w') as f:
        for line in sorted(results):
            f.write(line + '\n')

def main():
    parser = argparse.ArgumentParser(description="Analyze TSV file for pfam taxonomic distribution")
    parser.add_argument('--input', required=True, help="Input TSV file path")
    parser.add_argument('--output', required=True, help="Output TSV file path")
    
    args = parser.parse_args()
    
    analyze_tsv(args.input, args.output)
    print(f"Analysis complete. Results written to {args.output}")

if __name__ == "__main__":
    main()
