import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from tqdm import tqdm

def extract_protein_id(filename):
    """Extract protein ID from filename - first phrase before '__'"""
    return filename.split('__')[0]

def get_aa_chain_length(structure):
    """Count number of residues in the first chain of the structure"""
    for model in structure:
        for chain in model:
            return len(list(chain))
    return 0

def get_ca_coords(structure):
    """Extract CA coordinates from structure"""
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
    return np.array(coords)

def max_distance(coords):
    """Calculate maximum distance between CA atoms"""
    if len(coords) < 2:
        return 0.0
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))
    return np.max(dist_matrix)

def read_group_files(group_dir):
    """Read group files from directory and map protein IDs to group names and colors"""
    protein_to_group = {}
    
    # Define specific colors for eBAP groups 1-5
    group_to_color = {
        'eBAP1': '#3F78C1',  # Blue
        'eBAP2': '#9163C1',  # Purple
        'eBAP3': '#44A043',  # Green
        'eBAP4': '#D37538',  # Orange
        'eBAP5': '#B33D90'   # Magenta/Pink
    }
    
    if not os.path.exists(group_dir):
        print(f"WARNING: Group directory '{group_dir}' does not exist!")
        return {}, {}
    
    print(f"\n=== GROUP FILE ANALYSIS ===")
    print(f"Reading group files from: {group_dir}")
    
    group_files = [f for f in os.listdir(group_dir) if os.path.isfile(os.path.join(group_dir, f))]
    print(f"Found {len(group_files)} files in group directory: {group_files}")
    
    for filename in group_files:
        filepath = os.path.join(group_dir, filename)
        group_name = os.path.splitext(filename)[0]  # Use filename without extension as group name
        
        # Use predefined color if it's an eBAP group, otherwise use gray
        color = group_to_color.get(group_name, 'lightgray')
        
        print(f"\nProcessing group file: {filename}")
        print(f"Group name: {group_name}, Color: {color}")
        
        try:
            with open(filepath, 'r') as f:
                protein_ids = f.read().splitlines()
                valid_ids = [pid.strip() for pid in protein_ids if pid.strip()]
                print(f"Found {len(valid_ids)} protein IDs in {filename}:")
                for pid in valid_ids:
                    print(f"  - {pid}")
                    protein_to_group[pid] = group_name
        except Exception as e:
            print(f"ERROR reading {filename}: {e}")
    
    print(f"\nTotal protein IDs mapped to groups: {len(protein_to_group)}")
    print("="*50)
    
    return protein_to_group, group_to_color

def extract_fiber_data(pdb_path):
    """Extract both fiber length and AA chain length from PDB file"""
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('fiber', pdb_path)
        coords = get_ca_coords(structure)
        fiber_length = max_distance(coords)
        aa_chain_length = get_aa_chain_length(structure)
        return fiber_length, aa_chain_length
    except Exception as e:
        return np.nan, np.nan

def remove_outliers_iqr(df, column):
    """Remove outliers using IQR method"""
    Q1 = df[column].quantile(0.25)
    Q3 = df[column].quantile(0.75)
    IQR = Q3 - Q1
    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    df_clean = df[(df[column] >= lower) & (df[column] <= upper)]
    print(f"Removed {len(df) - len(df_clean)} outliers from '{column}' using IQR method.")
    return df_clean

def analyze_protein_id_matching(df, protein_to_group):
    """Analyze and report protein ID matching"""
    print(f"\n=== PROTEIN ID MATCHING ANALYSIS ===")
    
    unique_protein_ids = set(df['protein_id'].unique())
    group_protein_ids = set(protein_to_group.keys())
    
    print(f"Unique protein IDs extracted from PDB files ({len(unique_protein_ids)}):")
    for pid in sorted(unique_protein_ids):
        print(f"  - {pid}")
    
    print(f"\nProtein IDs from group files ({len(group_protein_ids)}):")
    for pid in sorted(group_protein_ids):
        print(f"  - {pid}")
    
    matched_ids = unique_protein_ids.intersection(group_protein_ids)
    unmatched_pdb_ids = unique_protein_ids - group_protein_ids
    unmatched_group_ids = group_protein_ids - unique_protein_ids
    
    print(f"\nMATCHED protein IDs ({len(matched_ids)}):")
    for pid in sorted(matched_ids):
        group = protein_to_group[pid]
        print(f"  - {pid} -> Group: {group}")
    
    print(f"\nUNMATCHED protein IDs from PDB files ({len(unmatched_pdb_ids)}):")
    for pid in sorted(unmatched_pdb_ids):
        print(f"  - {pid} (will be colored as 'Ungrouped')")
    
    print(f"\nUnused protein IDs from group files ({len(unmatched_group_ids)}):")
    for pid in sorted(unmatched_group_ids):
        print(f"  - {pid}")
    
    print("="*50)
    
    return matched_ids, unmatched_pdb_ids

def create_scatter_plot(df, output_file, protein_to_group=None, group_to_color=None):
    """Create scatter plot of AA chain length vs fiber length"""
    plt.figure(figsize=(12, 8))
    
    if protein_to_group and group_to_color:
        # Add group column to dataframe
        df['group'] = df['protein_id'].map(protein_to_group).fillna('Ungrouped')
        
        # Define the ordered sequence for eBAP groups
        ordered_groups = ['eBAP1', 'eBAP2', 'eBAP3', 'eBAP4', 'eBAP5']
        
        # Plot eBAP groups in order first
        for group_name in ordered_groups:
            if group_name in df['group'].values:
                group_data = df[df['group'] == group_name]
                color = group_to_color.get(group_name, 'lightgray')
                
                plt.scatter(group_data['aa_chain_length'], group_data['fiber_length'], 
                           label=f'{group_name} (n={len(group_data)})', 
                           color=color, alpha=0.7, s=60)
                
                print(f"Plotted group '{group_name}': {len(group_data)} points, color: {color}")
        
        # Plot any remaining groups (like 'Ungrouped' or other groups)
        remaining_groups = [g for g in df['group'].unique() if g not in ordered_groups]
        for group_name in remaining_groups:
            group_data = df[df['group'] == group_name]
            color = group_to_color.get(group_name, 'lightgray')
            
            plt.scatter(group_data['aa_chain_length'], group_data['fiber_length'], 
                       label=f'{group_name} (n={len(group_data)})', 
                       color=color, alpha=0.7, s=60)
            
            print(f"Plotted group '{group_name}': {len(group_data)} points, color: {color}")
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        # Single color for all points
        plt.scatter(df['aa_chain_length'], df['fiber_length'], 
                   color='blue', alpha=0.7, s=60, label=f'All data (n={len(df)})')
        plt.legend()
    
    plt.xlabel('AA Chain Length')
    plt.ylabel('Fiber Length (?)')
    plt.title('Scatter Plot of AA Chain Length vs Fiber Length')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, format='svg', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Extract fiber data from PDBs and create scatter plot.')
    parser.add_argument('input_folder', help='Folder containing PDB files (recursively searched)')
    parser.add_argument('output_tsv', help='Output TSV file name')
    parser.add_argument('output_scatter', help='Output scatter plot SVG file name (e.g., scatter.svg)')
    parser.add_argument('--remove_outliers', type=int, default=0, help='Set to 1 to remove fiber length outliers using IQR method')
    parser.add_argument('--groups_dir', help='Optional directory containing group files with protein IDs for color coding')
    args = parser.parse_args()

    # Read group information if directory provided
    protein_to_group = {}
    group_to_color = {}
    if args.groups_dir:
        protein_to_group, group_to_color = read_group_files(args.groups_dir)

    # Find all PDB files
    pdb_files = []
    for root, dirs, files in os.walk(args.input_folder):
        for file in files:
            if file.lower().endswith('.pdb'):
                pdb_files.append(os.path.join(root, file))

    print(f"\nFound {len(pdb_files)} PDB files to process")

    # Process PDB files
    records = []
    with tqdm(total=len(pdb_files), desc="Processing PDBs", unit="file") as pbar:
        for pdb_path in pdb_files:
            filename = os.path.basename(pdb_path)
            protein_id = extract_protein_id(filename)
            fiber_length, aa_chain_length = extract_fiber_data(pdb_path)
            
            records.append({
                'protein_id': protein_id,
                'pdb_file': filename, 
                'aa_chain_length': aa_chain_length,
                'fiber_length': fiber_length
            })
            
            pbar.set_postfix({
                "protein": protein_id,
                "AA_len": f"{aa_chain_length}" if not np.isnan(aa_chain_length) else "NaN",
                "fiber_len": f"{fiber_length:.2f}?" if not np.isnan(fiber_length) else "NaN"
            })
            pbar.update(1)

    # Create DataFrame and remove NaN values
    df = pd.DataFrame(records).dropna()
    
    if len(df) == 0:
        print("No valid data found. Please check your PDB files.")
        return

    # Analyze protein ID matching if groups are provided
    if protein_to_group:
        matched_ids, unmatched_ids = analyze_protein_id_matching(df, protein_to_group)

    # Remove outliers if requested
    if args.remove_outliers:
        df = remove_outliers_iqr(df, 'fiber_length')

    # Add group information to the dataframe for TSV output
    if protein_to_group:
        df['group'] = df['protein_id'].map(protein_to_group).fillna('Ungrouped')

    # Save TSV
    df.to_csv(args.output_tsv, sep='\t', index=False)
    
    # Create scatter plot
    create_scatter_plot(df, args.output_scatter, 
                       protein_to_group if protein_to_group else None,
                       group_to_color if group_to_color else None)
    
    print(f"\n=== PROCESSING COMPLETE ===")
    print(f"TSV saved to: {args.output_tsv}")
    print(f"Scatter plot saved to: {args.output_scatter}")
    print(f"Processed {len(df)} valid structures")
    
    if protein_to_group:
        print(f"\nGroup summary:")
        group_counts = df['group'].value_counts()
        for group, count in group_counts.items():
            color = group_to_color.get(group, 'lightgray')
            print(f"  - {group}: {count} structures (color: {color})")

if __name__ == '__main__':
    main()
