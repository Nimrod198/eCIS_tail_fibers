import csv
from Bio import Entrez
import time
import sys
from tqdm.auto import tqdm

def read_taxids_from_file(file_path):
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]

def get_taxonomy_info(taxids):
    taxonomy_info = {}
    taxids_not_found = []
    batch_size = 100

    with tqdm(total=len(taxids), desc="Processing TaxIDs", unit="taxid") as pbar:
        for i in range(0, len(taxids), batch_size):
            batch = taxids[i:i+batch_size]
            try:
                handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
                records = Entrez.read(handle)
                handle.close()

                for record in records:
                    taxid = record["TaxId"]
                    lineage = record.get("LineageEx", [])
                    lineage_dict = {item["Rank"]: item["ScientificName"] for item in lineage}
                    lineage_dict["species"] = record.get("ScientificName", "N/A")
                    taxonomy_info[taxid] = lineage_dict

                pbar.update(len(batch))
                time.sleep(1)  # Wait between batches
            except Exception as e:
                taxids_not_found.extend(batch)
                pbar.update(len(batch))

    if taxids_not_found:
        print(f"\nCould not retrieve information for {len(taxids_not_found)} TaxIDs")

    return taxonomy_info, taxids_not_found

def main(input_file, output_file):
    Entrez.email = "your_email@example.com"  # Replace with your email
    taxids = read_taxids_from_file(input_file)

    print(f"Read {len(taxids)} TaxIDs from {input_file}")

    taxonomy_info, taxids_not_found = get_taxonomy_info(taxids)

    ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

    print(f"Writing results to {output_file}")
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(["TaxID"] + ranks)
        
        with tqdm(total=len(taxonomy_info), desc="Writing results", unit="row") as pbar:
            for taxid, lineage in taxonomy_info.items():
                row = [taxid] + [lineage.get(rank, "N/A") for rank in ranks]
                writer.writerow(row)
                pbar.update(1)

    print(f"Wrote {len(taxonomy_info)} rows to {output_file}")

    if taxids_not_found:
        print(f"{len(taxids_not_found)} TaxIDs could not be retrieved from NCBI")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        main(input_file, output_file)
