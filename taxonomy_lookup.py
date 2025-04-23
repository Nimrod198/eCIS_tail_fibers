import re
import csv
import logging
from Bio import SeqIO
from Bio import Entrez
import time

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def extract_taxids_and_proteins(fasta_file):
    taxid_pattern = re.compile(r"taxID:(\d+)")
    protein_taxid_map = {}
    proteins_without_taxid = []
    
    with open(fasta_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            header = record.description
            protein_id = record.id.split("|")[0]  # Get the first part of the ID before the |
            
            match = taxid_pattern.search(header)
            if match:
                taxid = match.group(1)
                protein_taxid_map[protein_id] = taxid
                logging.info(f"Found TaxID {taxid} for protein {protein_id}")
            else:
                proteins_without_taxid.append(protein_id)
                logging.warning(f"No TaxID found for protein {protein_id}")
    
    logging.info(f"Processed {len(protein_taxid_map) + len(proteins_without_taxid)} proteins")
    logging.info(f"Found TaxIDs for {len(protein_taxid_map)} proteins")
    logging.info(f"No TaxID found for {len(proteins_without_taxid)} proteins")
    
    return protein_taxid_map, proteins_without_taxid

def get_taxonomy_info(taxids):
    taxonomy_info = {}
    taxids_not_found = []
    batch_size = 100
    
    for i in range(0, len(taxids), batch_size):
        batch = list(taxids)[i:i+batch_size]
        try:
            handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            
            for record in records:
                taxid = record["TaxId"]
                lineage = record.get("Lineage", "N/A").split("; ")
                
                first_rank = lineage[2] if len(lineage) > 2 else "N/A"
                last_rank = lineage[-1] if lineage else "N/A"
                
                taxonomy_info[taxid] = {
                    "first_rank": first_rank,
                    "last_rank": last_rank
                }
                logging.info(f"Retrieved taxonomy info for TaxID {taxid}")
            
            logging.info(f"Processed batch of {len(batch)} TaxIDs")
            time.sleep(1)  # Wait between batches
        except Exception as e:
            logging.error(f"Error fetching batch: {str(e)}")
            taxids_not_found.extend(batch)
    
    if taxids_not_found:
        logging.warning(f"Could not retrieve information for {len(taxids_not_found)} TaxIDs")
    
    return taxonomy_info, taxids_not_found

def main(fasta_file, output_file):
    protein_taxid_map, proteins_without_taxid = extract_taxids_and_proteins(fasta_file)
    unique_taxids = set(protein_taxid_map.values())
    
    logging.info(f"Extracted {len(protein_taxid_map)} protein IDs with TaxIDs and {len(proteins_without_taxid)} without TaxIDs")
    logging.info(f"Found {len(unique_taxids)} unique TaxIDs")
    
    taxonomy_info, taxids_not_found = get_taxonomy_info(unique_taxids)
    
    logging.info(f"Writing results to {output_file}")
    rows_written = 0
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Protein_ID", "TaxID", "First_Rank", "Last_Rank"])
        for protein_id, taxid in protein_taxid_map.items():
            if taxid in taxonomy_info:
                first_rank = taxonomy_info[taxid]["first_rank"]
                last_rank = taxonomy_info[taxid]["last_rank"]
                writer.writerow([protein_id, taxid, first_rank, last_rank])
                rows_written += 1
            else:
                logging.warning(f"No taxonomy information found for TaxID {taxid} (Protein ID: {protein_id})")
    
    logging.info(f"Wrote {rows_written} rows to {output_file}")
    
    if proteins_without_taxid:
        logging.warning(f"{len(proteins_without_taxid)} proteins had no TaxID and were not included in the output")
    
    if taxids_not_found:
        logging.warning(f"{len(taxids_not_found)} TaxIDs could not be retrieved from NCBI")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fasta_file> <output_file>")
    else:
        fasta_file = sys.argv[1]
        output_file = sys.argv[2]
        main(fasta_file, output_file)