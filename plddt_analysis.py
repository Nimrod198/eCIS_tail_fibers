import os
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

def extract_avg_plddt(pdb_file):
    """Extract average pLDDT score from a PDB file"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('', pdb_file)
    plddt_scores = [atom.get_bfactor() 
                    for atom in structure.get_atoms()]
    avg_plddt = sum(plddt_scores)/len(plddt_scores) if plddt_scores else 0
    file_id = os.path.splitext(os.path.basename(pdb_file))[0]
    return file_id, avg_plddt

def analyze_pdb_dir(directory):
    """Process directory of PDB files and generate outputs"""
    results = []
    for filename in os.listdir(directory):
        if filename.endswith('.pdb'):
            file_path = os.path.join(directory, filename)
            file_id, avg_plddt = extract_avg_plddt(file_path)
            results.append({'ID': file_id, 'Average_pLDDT': avg_plddt})
    
    # Create and save data table
    df = pd.DataFrame(results)
    csv_path = os.path.join(directory, 'plddt_summary.csv')
    df.to_csv(csv_path, index=False)
    
    # Generate plots
    plt.figure(figsize=(10, 6))
    plt.boxplot(df['Average_pLDDT'], vert=True)
    plt.title('Average pLDDT Distribution (Boxplot)')
    plt.ylabel('pLDDT Score')
    plt.savefig(os.path.join(directory, 'plddt_boxplot.png'))
    plt.close()
    
    plt.figure(figsize=(10, 6))
    plt.violinplot(df['Average_pLDDT'], showmeans=True)
    plt.title('Average pLDDT Distribution (Violin Plot)')
    plt.ylabel('pLDDT Score')
    plt.savefig(os.path.join(directory, 'plddt_violin.png'))
    plt.close()
    
    return csv_path

# Usage:
# analyze_pdb_dir('/path/to/pdb_directory')
