import os
import csv
import argparse
import logging
from datetime import datetime
from tqdm import tqdm

def setup_logging(log_file):
    """Set up logging with both console and file handlers."""
    logger = logging.getLogger('gff_search_logger')
    logger.setLevel(logging.DEBUG)  # Set to DEBUG to capture all levels

    # Console Handler (INFO and above)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_formatter)

    # File Handler (DEBUG and above)
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')
    file_handler.setFormatter(file_formatter)

    # Add handlers to logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger

def read_tsv_file(tsv_file, logger):
    """Read the input TSV file and extract genome and gene IDs along with Pfam data."""
    pfam_data = []
    logger.info(f"Reading input TSV file: {tsv_file}")
    
    try:
        with open(tsv_file, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                pfam_data.append({
                    'genomeID': row['genomeID'],
                    'geneID': row['geneID'],
                    'pfamID': row['pfamID'],
                    'pfamDescription': row['pfamDescription']
                })
        logger.info(f"Total entries read from TSV: {len(pfam_data)} genes")
    except Exception as e:
        logger.error(f"Error reading TSV file: {e}")
    
    return pfam_data

def search_gff_for_genes(gff_file, gene_ids, logger):
    """Search a GFF file for entries matching the given gene IDs."""
    gff_data = []
    logger.debug(f"Searching for gene IDs: {gene_ids}")
    
    try:
        with open(gff_file, 'r') as file:
            for line_num, line in enumerate(file, 1):
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                    seqid = parts[0]
                    attributes = parts[8]
                    try:
                        attr_dict = dict(item.split('=') for item in attributes.split(';') if '=' in item)
                        gff_id = attr_dict.get('ID', '')
                        
                        for gene_id in gene_ids:
                            if '_' in gene_id:
                                # Complex ID handling
                                gene_prefix, gene_suffix = gene_id.split('_', 1)
                                if seqid == gene_prefix and '_' in gff_id:
                                    gff_suffix = gff_id.split('_', 1)[1]
                                    if gene_suffix == gff_suffix:
                                        logger.debug(f"Match found for complex ID: {gene_id}")
                                        gff_data.append(parts + [gene_id])
                                        break
                            else:
                                # Standard ID handling
                                if gff_id == gene_id:
                                    logger.debug(f"Match found for standard ID: {gene_id}")
                                    gff_data.append(parts + [gene_id])
                                    break
                    except ValueError as e:
                        logger.warning(f"Error parsing attributes in line {line_num}: {line.strip()}")
                        logger.warning(f"Exception: {e}")
    except IOError as e:
        logger.error(f"Error reading GFF file {gff_file}: {e}")
    
    logger.debug(f"Total matching entries found: {len(gff_data)}")
    return gff_data

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Search GFF files for genes and merge with Pfam data.')
    parser.add_argument('--input-tsv', required=True, help='Input TSV file containing genome and gene IDs with Pfams.')
    parser.add_argument('--gff-directory', required=True, help='Directory containing GFF files.')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file to save results.')
    parser.add_argument('--log-file', 
                        default=f'gff_search_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log', 
                        help='Log file to record processing details.')

    args = parser.parse_args()

    # Set up logging
    logger = setup_logging(args.log_file)

    # Print and log initial processing information
    logger.info("=" * 50)
    logger.info("GFF Search and Merge Script Started")
    logger.info("=" * 50)
    logger.info(f"Input TSV: {args.input_tsv}")
    logger.info(f"GFF Directory: {args.gff_directory}")
    logger.info(f"Output TSV: {args.output_tsv}")  # Corrected from output_tsV to output_tsv
    logger.info(f"Log File: {args.log_file}")
    logger.info("-" * 50)

    # Read the input TSV file
    pfam_data = read_tsv_file(args.input_tsv, logger)

    # Prepare tracking variables
    missing_gff_files = []
    processed_entries = 0
    matching_entries = 0
    genome_match_counts = {}
    
    # Open output file and write header
    with open(args.output_tsv, 'w', newline='') as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        
        # Write header
        writer.writerow(['genomeID', 'geneID', 'pfamID', 'pfamDescription', 
                         'seqid', 'source', 'type', 'start', 'end', 
                         'score', 'strand', 'phase', 'attributes', 'gene_id'])

        # Process each genome with a progress bar
        for entry in tqdm(pfam_data, desc="Processing Genes", unit="gene"):
            genome_id = entry['genomeID']
            gene_id = entry['geneID']
            gff_file_path = os.path.join(args.gff_directory, f"{genome_id}.gff")
            
            if os.path.exists(gff_file_path):
                # Search for matching genes in GFF
                gff_entries = search_gff_for_genes(gff_file_path, {gene_id}, logger)
                
                if gff_entries:
                    # Log detailed information about matching entries
                    logger.info(f"Genome {genome_id}, Gene {gene_id}: Found {len(gff_entries)} matching entries")
                    
                    # Track genome-specific match counts
                    genome_match_counts[genome_id] = genome_match_counts.get(genome_id, 0) + len(gff_entries)
                    
                    # Write matching entries
                    for gff_entry in gff_entries:
                        combined_row = [genome_id, gene_id, entry['pfamID'], entry['pfamDescription']] + gff_entry
                        writer.writerow(combined_row)
                        matching_entries += 1

                    output_file.flush()
                
                processed_entries += 1
            else:
                # Log missing GFF files
                logger.warning(f"GFF file not found for Genome {genome_id}: {gff_file_path}")
                missing_gff_files.append(gff_file_path)

    # Final summary logging
    logger.info("\n" + "=" * 50)
    logger.info("Processing Summary")
    logger.info("=" * 50)
    logger.info(f"Total genes processed: {processed_entries}")
    logger.info(f"Total matching entries written: {matching_entries}")

    # Log genome-specific match counts
    logger.info("\nGenome-specific Match Counts:")
    for genome, count in sorted(genome_match_counts.items(), key=lambda x: x[1], reverse=True):
        logger.info(f"  {genome}: {count} matching entries")
    
    # Log missing GFF files
    if missing_gff_files:
        logger.warning("\nMissing GFF Files:")
        for file_path in missing_gff_files:
            logger.warning(f"  {file_path}")
        logger.warning(f"Total missing GFF files: {len(missing_gff_files)}")

    logger.info("\nGFF Search and Merge Script Completed")

if __name__ == '__main__':
     main()