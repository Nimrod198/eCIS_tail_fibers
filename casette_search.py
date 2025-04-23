import os
import argparse
from tqdm import tqdm

def parse_pfam_file(pfam_file):
    """Parse a Pfam file and return gene and genome details."""
    gene_details = []
    with open(pfam_file, 'r') as file:
        header = file.readline()  # Read and skip the header line
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 12:  # Ensure there are enough columns
                pfam_id_full = parts[0].rstrip('.')  # Remove any trailing dot
                pfam_id = pfam_id_full.split('.')[0]  # Take only the ID part (before the dot)
                gene_id = parts[3]  # Query_ID
                pfam_description = parts[-1]  # Description
                gene_details.append((gene_id, pfam_id, pfam_description))
    return gene_details

def read_genome_list(genome_list_file):
    """Read genome IDs from a file."""
    with open(genome_list_file, 'r') as file:
        genomes = [line.strip() for line in file if line.strip()]
    return genomes

def search_pfam_in_genomes(pfam_ids, genome_ids, directory):
    """Search for specified Pfam IDs in the given genome files."""
    results = []
    total_genomes = len(genome_ids)
    total_pfams_retrieved = 0
    
    for index, genome_id in enumerate(tqdm(genome_ids, desc="Processing Genomes", unit="genome")):
        pfam_file_path = os.path.join(directory, f"{genome_id}.pfam")
        if os.path.exists(pfam_file_path):
            gene_details = parse_pfam_file(pfam_file_path)
            for gene_id, pfam_id, pfam_description in gene_details:
                if pfam_id in pfam_ids:
                    results.append((genome_id, gene_id, pfam_id, pfam_description))
                    total_pfams_retrieved += 1
        else:
            print(f"Warning: {pfam_file_path} does not exist.")
    
    print(f"\nSummary:")
    print(f"Total Genomes Processed: {total_genomes}")
    print(f"Total Pfams Retrieved: {total_pfams_retrieved}")

    return results

def write_results_to_file(results, output_file):
    """Write the results to a specified output file in TSV format."""
    if results:  # Check if there are results to write
        with open(output_file, 'w') as file:
            file.write("genomeID\tgeneID\tpfamID\tpfamDescription\n")
            for genome_id, gene_id, pfam_id, pfam_description in results:
                file.write(f"{genome_id}\t{gene_id}\t{pfam_id}\t{pfam_description}\n")
        print(f"Results written to {output_file}")
    else:
        print("No results to write.")

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Search Pfam files for specified Pfams in given genomes.')
    parser.add_argument('--pfams', required=True, help='Comma-separated list of Pfam IDs to search for.')
    parser.add_argument('--genomes', required=True, help='Path to the file containing genome IDs to search in.')
    parser.add_argument('--directory', required=True, help='Directory containing Pfam files.')
    parser.add_argument('--output', required=True, help='Output file to save results.')

    args = parser.parse_args()

    # Split the Pfams into a list
    pfam_ids = args.pfams.split(',')
    
    # Read the genome list from the specified file
    genome_ids = read_genome_list(args.genomes)

    # Perform the search
    results = search_pfam_in_genomes(pfam_ids, genome_ids, args.directory)

    # Write the results to an output file in TSV format
    write_results_to_file(results, args.output)

if __name__ == '__main__':
    main()