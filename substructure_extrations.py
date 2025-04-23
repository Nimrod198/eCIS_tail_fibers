import os
import pandas as pd
import argparse
from Bio import PDB
from tqdm import tqdm  # Import tqdm for progress bar

def extract_region(structure, chain_id, start, end):
    """Extract a specific region from a given chain based on start and end coordinates."""
    chain = structure[0][chain_id]
    new_chain = PDB.Chain.Chain(chain_id)
    for residue in chain:
        if start <= residue.id[1] <= end:
            new_chain.add(residue.copy())
    return new_chain

def process_structure(input_file, accession, start, end, output_dir):
    """Process the structure file to extract the specified region."""
    output_file = os.path.join(output_dir, f"{accession}_{start}_{end}.pdb")
    
    # Check if the file already exists
    if os.path.exists(output_file):
        return

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(accession, input_file)
    
    new_structure = PDB.Structure.Structure(accession)
    model = PDB.Model.Model(0)
    new_structure.add(model)
    
    for chain in structure[0]:
        extracted_chain = extract_region(structure, chain.id, start, end)
        if len(extracted_chain) > 0:
            model.add(extracted_chain)
    
    # Save the extracted region as a PDB file
    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(output_file)

def main():
    parser = argparse.ArgumentParser(description="Extract specific regions from protein structures based on DSS score.")
    parser.add_argument("input_tsv", help="Path to the input TSV file.")
    parser.add_argument("structure_dir", help="Directory containing the structure files.")
    parser.add_argument("output_dir", help="Directory to save the extracted structures.")
    parser.add_argument("--min-dss", type=float, default=0, help="Minimum DSS score required for extraction")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Read the input TSV file
    df = pd.read_csv(args.input_tsv, sep="\t")

    # Ensure columns are correctly named
    df.columns = ['Gene', 'Domain', 'HitCount', 'Coverage', 'AvgBitScore', 'DSS']

    # Extract start and end coordinates from Domain column
    df[['start_coordinate', 'end_coordinate']] = df['Domain'].str.split('-', expand=True).astype(int)

    # Filter rows based on minimum DSS score
    filtered_df = df[df['DSS'] >= args.min_dss]

    # Initialize progress bar
    with tqdm(total=len(filtered_df), desc="Processing structures") as pbar:
        for _, row in filtered_df.iterrows():
            # Extract accession from the Gene column
            accession = row["Gene"].split("__")[0]
            
            start = row["start_coordinate"]   # Start coordinate
            end = row["end_coordinate"]       # End coordinate
            
            # Find structure file in the specified directory
            matching_files = [f for f in os.listdir(args.structure_dir) if accession in f and f.endswith(".pdb")]
            
            if matching_files:
                input_file = os.path.join(args.structure_dir, matching_files[0])  # Select the first matching file
                try:
                    process_structure(input_file, accession, start, end, args.output_dir)
                    pbar.set_postfix({'status': f"Processed {accession}"})  # Update status in the progress bar
                except Exception as e:
                    print(f"Error processing {accession}: {str(e)}")
            else:
                print(f"Structure file not found for {accession}")
                pbar.set_postfix({'status': f"No file for {accession}"})  # Update status in the progress bar

            # Update progress bar
            pbar.update(1)

if __name__ == "__main__":
    main()
