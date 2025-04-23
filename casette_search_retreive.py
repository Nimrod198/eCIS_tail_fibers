import os
import csv
import argparse
import logging
from datetime import datetime
from tqdm import tqdm

def setup_logging(log_file):
    """Set up logging with both console and file handlers."""
    logger = logging.getLogger('gff_search_logger')
    logger.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_formatter)

    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')
    file_handler.setFormatter(file_formatter)

    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    return logger

def read_tsv_file(tsv_file, logger):
    """Read the input TSV file with cluster ID, genome number, and gene ID."""
    cluster_data = []
    try:
        with open(tsv_file, 'r') as file:
            for line in file:
                cluster_id, genome_number, gene_id = line.strip().split('\t')
                cluster_data.append({
                    'clusterID': cluster_id,
                    'genomeNumber': genome_number,
                    'geneID': gene_id
                })
        logger.info(f"Total entries read from TSV: {len(cluster_data)} genes")
    except Exception as e:
        logger.error(f"Error reading TSV file: {e}")
    return cluster_data

def search_gff_for_genes_hybrid(gff_file, gene_id, logger):
    """
    Hybrid search function that can handle both simple and complex gene IDs.
    - For simple numeric IDs: uses exact matching (original script behavior)
    - For complex IDs with underscores: uses the splitting logic to match parts
    """
    gff_data = []
    logger.debug(f"Searching for gene ID: {gene_id}")
    
    # Determine if this is a complex ID that needs splitting
    is_complex_with_underscore = '_' in gene_id
    
    if is_complex_with_underscore:
        # Complex ID case - split the gene ID for matching parts
        supplied_seqid, supplied_last_number = gene_id.rsplit('_', 1)
        logger.debug(f"Using complex ID matching with split: {supplied_seqid} and {supplied_last_number}")
        matching_mode = "complex"
    else:
        # Simple ID case - use original exact matching behavior
        logger.debug(f"Using simple exact matching for ID: {gene_id}")
        matching_mode = "simple"
    
    try:
        with open(gff_file, 'r') as file:
            for line_num, line in enumerate(file, 1):
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                gff_seqid = parts[0]
                attributes = parts[8]
                
                # Parse attributes
                try:
                    attr_dict = dict(item.split('=') for item in attributes.split(';') if '=' in item)
                    gff_id = attr_dict.get('ID', '')
                except ValueError as e:
                    logger.warning(f"Error parsing attributes in line {line_num}: {line.strip()}")
                    logger.warning(f"Exception: {e}")
                    continue
                
                # Apply appropriate matching strategy based on ID type
                if matching_mode == "complex":
                    # For complex IDs, match both parts (seqid and suffix)
                    if '_' in gff_id:
                        _, gff_last_number = gff_id.rsplit('_', 1)
                        if gff_seqid == supplied_seqid and gff_last_number == supplied_last_number:
                            logger.debug(f"Complex match found in line {line_num}: seqid={gff_seqid}, number={gff_last_number}")
                            gff_data.append(parts + [gene_id])
                else:
                    # For simple IDs, do exact matching on ID attribute (original behavior)
                    if gff_id == gene_id:
                        logger.debug(f"Simple match found in line {line_num}: ID={gff_id}")
                        gff_data.append(parts + [gene_id])
    except IOError as e:
        logger.error(f"Error reading GFF file {gff_file}: {e}")
    
    logger.debug(f"Total matching entries found for gene ID {gene_id}: {len(gff_data)}")
    return gff_data

def main():
    parser = argparse.ArgumentParser(description='Search GFF files for genes and merge with cluster data.')
    parser.add_argument('--input-tsv', required=True, help='Input TSV file containing cluster ID, genome number, and gene ID.')
    parser.add_argument('--gff-directory', required=True, help='Directory containing GFF files.')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file to save results.')
    parser.add_argument('--log-file', default=f'gff_search_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log', help='Log file to record processing details.')
    args = parser.parse_args()

    logger = setup_logging(args.log_file)

    logger.info("=" * 50)
    logger.info("GFF Search and Merge Script Started")
    logger.info("=" * 50)
    logger.info(f"Input TSV: {args.input_tsv}")
    logger.info(f"GFF Directory: {args.gff_directory}")
    logger.info(f"Output TSV: {args.output_tsv}")
    logger.info(f"Log File: {args.log_file}")
    logger.info("-" * 50)

    cluster_data = read_tsv_file(args.input_tsv, logger)

    missing_gff_files = []
    processed_entries = 0
    matching_entries = 0
    genome_match_counts = {}

    with open(args.output_tsv, 'w', newline='') as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerow(['clusterID', 'genomeNumber', 'geneID', 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

        for entry in tqdm(cluster_data, desc="Processing Genes", unit="gene"):
            cluster_id = entry['clusterID']
            genome_number = entry['genomeNumber']
            gene_id = entry['geneID']
            gff_file_path = os.path.join(args.gff_directory, f"{genome_number}.gff")

            if os.path.exists(gff_file_path):
                gff_entries = search_gff_for_genes_hybrid(gff_file_path, gene_id, logger)
                if gff_entries:
                    logger.info(f"Genome {genome_number}, Gene {gene_id}: Found {len(gff_entries)} matching entries")
                    genome_match_counts[genome_number] = genome_match_counts.get(genome_number, 0) + len(gff_entries)
                    for gff_entry in gff_entries:
                        combined_row = [cluster_id, genome_number, gene_id] + gff_entry[:-1]
                        writer.writerow(combined_row)
                        matching_entries += 1
                    output_file.flush()
                processed_entries += 1
            else:
                logger.warning(f"GFF file not found for Genome {genome_number}: {gff_file_path}")
                missing_gff_files.append(gff_file_path)

    logger.info("\n" + "=" * 50)
    logger.info("Processing Summary")
    logger.info("=" * 50)
    logger.info(f"Total genes processed: {processed_entries}")
    logger.info(f"Total matching entries written: {matching_entries}")

    logger.info("\nGenome-specific Match Counts:")
    for genome, count in sorted(genome_match_counts.items(), key=lambda x: x[1], reverse=True):
        logger.info(f" {genome}: {count} matching entries")

    if missing_gff_files:
        logger.warning("\nMissing GFF Files:")
        for file_path in missing_gff_files:
            logger.warning(f" {file_path}")
        logger.warning(f"Total missing GFF files: {len(missing_gff_files)}")

    logger.info("\nGFF Search and Merge Script Completed")

if __name__ == '__main__':
    main()
