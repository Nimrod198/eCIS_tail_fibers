import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def calculate_bulkiness_factor(df):
    """
    Calculate bulkiness factor = aa_chain_length / fiber_length
    Higher values indicate more bulky (compact) structures
    Lower values indicate more elongated structures
    """
    # Avoid division by zero
    df['bulkiness_factor'] = np.where(df['fiber_length'] > 0, 
                                     df['aa_chain_length'] / df['fiber_length'], 
                                     np.nan)
    return df

def create_bulkiness_boxplot(df, output_file):
    """Create boxplot showing bulkiness factor distribution by eBAP group"""
    
    # Filter for eBAP groups 1-5
    valid_groups = ['eBAP1', 'eBAP2', 'eBAP3', 'eBAP4', 'eBAP5']
    df_filtered = df[df['group'].isin(valid_groups)].copy()
    
    if df_filtered.empty:
        print("No data found for eBAP groups 1-5")
        return
    
    # Remove NaN values
    df_filtered = df_filtered.dropna(subset=['bulkiness_factor'])
    
    # Define specific colors for each eBAP group
    group_colors = {
        'eBAP1': '#3F78C1',  # Blue
        'eBAP2': '#9163C1',  # Purple
        'eBAP3': '#44A043',  # Green
        'eBAP4': '#D37538',  # Orange
        'eBAP5': '#B33D90'   # Magenta/Pink
    }
    
    # Prepare data for boxplot
    boxplot_data = []
    labels = []
    colors = []
    
    print("\n=== BULKINESS ANALYSIS SUMMARY ===")
    print("Bulkiness Factor = AA Chain Length / Fiber Length (?)")
    print("Higher values = More compact/bulky structures")
    print("Lower values = More elongated structures\n")
    
    for group in valid_groups:
        group_data = df_filtered[df_filtered['group'] == group]['bulkiness_factor']
        
        if len(group_data) > 0:
            boxplot_data.append(group_data.values)
            labels.append(f'{group}\n(n={len(group_data)})')
            colors.append(group_colors[group])
            
            # Print statistics
            print(f"{group} (n={len(group_data)}):")
            print(f"  Mean: {group_data.mean():.3f}")
            print(f"  Median: {group_data.median():.3f}")
            print(f"  Std: {group_data.std():.3f}")
            print(f"  Range: {group_data.min():.3f} - {group_data.max():.3f}")
            print()
    
    # Create the boxplot
    fig, ax = plt.subplots(figsize=(7, 5))
    
    # Create boxplot
    bp = ax.boxplot(boxplot_data, labels=labels, patch_artist=True, 
                    notch=True, showmeans=True)
    
    # Apply specific colors to boxes
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Customize box plot elements
    for element in ['whiskers', 'fliers', 'medians', 'caps']:
        plt.setp(bp[element], color='black')
    
    # Customize means (orange diamonds)
    plt.setp(bp['means'], marker='D', markerfacecolor='orange', 
             markeredgecolor='black', markersize=6)
    
    # Customize the plot
    ax.set_xlabel('eBAP Group', fontsize=12, fontweight='bold')
    ax.set_ylabel('(AA Chain Length / Fiber Length)', fontsize=12, fontweight='bold')
    ax.set_title('Protein Bulkiness Factor Distribution by eBAP Group', fontsize=14, fontweight='bold')
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
    
    # Add annotation
    ax.text(0.02, 0.98, 'Higher values = More compact structures\nLower values = More elongated structures', 
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))
    
    
    # Customize appearance
    plt.xticks(rotation=0)
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Boxplot saved to: {output_file}")
    
    # Show the plot
    plt.show()
    plt.close()

def analyze_group_differences(df):
    """Perform statistical analysis of differences between groups"""
    valid_groups = ['eBAP1', 'eBAP2', 'eBAP3', 'eBAP4', 'eBAP5']
    df_filtered = df[df['group'].isin(valid_groups)].dropna(subset=['bulkiness_factor'])
    
    print("\n=== GROUP COMPARISON ===")
    
    # Calculate group means
    group_means = df_filtered.groupby('group')['bulkiness_factor'].agg(['mean', 'std', 'count'])
    print("Group Statistics:")
    print(group_means.round(3))
    
    # Find most and least bulky groups
    mean_bulkiness = df_filtered.groupby('group')['bulkiness_factor'].mean()
    most_bulky = mean_bulkiness.idxmax()
    least_bulky = mean_bulkiness.idxmin()
    
    print(f"\nMost bulky group: {most_bulky} (mean bulkiness: {mean_bulkiness[most_bulky]:.3f})")
    print(f"Least bulky group: {least_bulky} (mean bulkiness: {mean_bulkiness[least_bulky]:.3f})")
    print(f"Difference: {mean_bulkiness[most_bulky] - mean_bulkiness[least_bulky]:.3f}x")
    
    # Print color coding for reference
    group_colors = {
        'eBAP1': '#3F78C1',  # Blue
        'eBAP2': '#9163C1',  # Purple  
        'eBAP3': '#44A043',  # Green
        'eBAP4': '#D37538',  # Orange
        'eBAP5': '#B33D90'   # Magenta/Pink
    }
    
    print(f"\nColor coding:")
    for group in valid_groups:
        if group in mean_bulkiness.index:
            print(f"  {group}: {group_colors[group]} (mean bulkiness: {mean_bulkiness[group]:.3f})")

def main():
    parser = argparse.ArgumentParser(description='Analyze protein bulkiness factors by eBAP group')
    parser.add_argument('input_tsv', help='Input TSV file with protein data')
    parser.add_argument('output_plot', help='Output boxplot file (e.g., bulkiness_boxplot.png)')
    parser.add_argument('--show_outliers', action='store_true', help='Show individual data points as outliers')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input_tsv):
        print(f"Error: Input file '{args.input_tsv}' not found!")
        return
    
    # Read the TSV data
    try:
        df = pd.read_csv(args.input_tsv, sep='\t')
        print(f"Loaded {len(df)} protein structures from {args.input_tsv}")
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        return
    
    # Check required columns
    required_columns = ['protein_id', 'aa_chain_length', 'fiber_length', 'group']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: Missing required columns: {missing_columns}")
        return
    
    # Calculate bulkiness factor
    df = calculate_bulkiness_factor(df)
    
    # Create boxplot
    create_bulkiness_boxplot(df, args.output_plot)
    
    # Perform statistical analysis
    analyze_group_differences(df)
    
    # Optional: Save processed data with bulkiness factor
    output_tsv = args.input_tsv.replace('.tsv', '_with_bulkiness.tsv')
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Enhanced data saved to: {output_tsv}")

if __name__ == '__main__':
    main()
