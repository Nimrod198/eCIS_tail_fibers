import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_data(data_string):
    print("Starting data parsing...")
    lines = data_string.strip().split('\n')
    data = []
    all_taxonomies = set()
    skipped_lines = 0
    
    for line in lines:
        try:
            parts = line.split(maxsplit=1)
            if len(parts) != 2:
                skipped_lines += 1
                continue
            
            family = parts[0]
            distributions = parts[1].strip('{}').split(';')
            dist_dict = {}
            
            for item in distributions:
                if item.strip():
                    percentage, taxonomy = item.split(',')
                    percentage = float(percentage.strip('%').strip('>').strip())
                    dist_dict[taxonomy] = percentage
                    all_taxonomies.add(taxonomy)
            
            if dist_dict:
                data.append([family] + [dist_dict.get(tax, 0) for tax in sorted(all_taxonomies)])
            else:
                skipped_lines += 1
        
        except Exception as e:
            print(f"Error parsing line: {line}")
            print(f"Error message: {str(e)}")
            skipped_lines += 1
    
    columns = ['Family'] + sorted(all_taxonomies)
    df = pd.DataFrame(data, columns=columns)
    
    print(f"Data parsing complete. {len(df)} families processed. {skipped_lines} lines skipped.")
    print(f"Taxonomies found: {', '.join(sorted(all_taxonomies))}")
    
    return df

def main(input_file):
    print(f"Reading input file: {input_file}")
    
    with open(input_file, 'r') as file:
        data_string = file.read()
    
    df = parse_data(data_string)
    
    if df.empty:
        print("No valid data to process. Exiting.")
        return
    
    print("\nGenerating heatmap...")
    
    # Sort the DataFrame by the highest share of taxonomy
    df_sorted = df.sort_values(by=df.columns[1:].tolist(), ascending=False)
    
    plt.figure(figsize=(20, 15))
    sns.heatmap(df_sorted.set_index('Family'), cmap='YlOrRd', annot=False, cbar_kws={'label': 'Percentage'})
    plt.title('Taxonomic Distribution Heatmap')
    plt.tight_layout()
    plt.savefig('heatmap.png', dpi=300)
    plt.close()
    print("Heatmap saved as heatmap.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    main(input_file)
