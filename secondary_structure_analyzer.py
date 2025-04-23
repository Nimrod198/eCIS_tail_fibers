import os
import sys
from Bio import PDB
from Bio.PDB.DSSP import DSSP
import csv
from tqdm import tqdm  # For progress bar

def analyze_pdb(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", pdb_file)
    except Exception as e:
        print(f"\nError parsing {os.path.basename(pdb_file)}: {str(e)}")
        return []

    model = structure[0]
    
    try:
        dssp = DSSP(model, pdb_file)
    except Exception as e:
        print(f"\nDSSP failed for {os.path.basename(pdb_file)}: {str(e)}")
        return []

    results = []
    current_sec = ''
    current_length = 0

    for key in dssp.keys():
        ss = dssp[key][2]
        if ss == 'H':  # Alpha
            if current_sec == 'A':
                current_length += 1
            else:
                if current_sec:
                    results.append((current_sec, current_length))
                current_sec = 'A'
                current_length = 1
        elif ss == 'E':  # Beta
            if current_sec == 'B':
                current_length += 1
            else:
                if current_sec:
                    results.append((current_sec, current_length))
                current_sec = 'B'
                current_length = 1
        else:  # Loop
            if current_sec == 'L':
                current_length += 1
            else:
                if current_sec:
                    results.append((current_sec, current_length))
                current_sec = 'L'
                current_length = 1

    # Append the last recorded structure
    if current_sec:
        results.append((current_sec, current_length))

    return results

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_dir> <output.tsv> [-v]")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    verbose = '-v' in sys.argv

    # Get list of PDB files first for accurate progress tracking
    pdb_files = [f for f in os.listdir(input_dir) if f.endswith(".pdb")]
    if not pdb_files:
        print("No PDB files found in directory!")
        sys.exit(1)

    with open(output_file, 'w', newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(['file_name', 'sec_structure'])

        # Initialize progress bar
        with tqdm(total=len(pdb_files), unit='file', 
                 desc="Processing PDB files") as pbar:
            for filename in pdb_files:
                pdb_path = os.path.join(input_dir, filename)
                
                if verbose:
                    pbar.set_postfix(file=filename[:20])
                
                try:
                    results = analyze_pdb(pdb_path)
                    for sec_type, length in results:
                        writer.writerow([filename, f'{length}aa({sec_type})'])
                except Exception as e:
                    print(f"\nError processing {filename}: {str(e)}")
                
                pbar.update(1)  # Update progress bar

    print(f"\nAnalysis complete! Results saved to {output_file}")

if __name__ == "__main__":
    main()
