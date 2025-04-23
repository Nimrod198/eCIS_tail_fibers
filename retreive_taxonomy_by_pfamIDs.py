import argparse
import requests
import json
import time
from collections import defaultdict

def fetch_pfam_data(pfam_id):
    pfam_id = pfam_id.split('.')[0]  # Remove version number
    url = f"https://www.ebi.ac.uk/interpro/api/taxonomy/uniprot/entry/pfam/{pfam_id}/"
    headers = {'Accept': 'application/json'}
    
    all_results = []
    while url:
        try:
            response = requests.get(url, headers=headers, timeout=10)
            response.raise_for_status()
            data = response.json()
            all_results.extend(data['results'])
            url = data.get('next')  # Get next page URL if it exists
            if url:
                time.sleep(1)  # Rate limiting between pages
        except requests.exceptions.RequestException as e:
            return f"Error: {str(e)}"
    
    return {'count': data.get('count', 0), 'results': all_results}

def is_kingdom_level(node):
    """Identify if a node represents a kingdom-level classification"""
    kingdom_names = {
        "Bacteria": "Bacteria",
        "Archaea": "Archaea",
        "Eukaryota": "Eukaryota",
        "Viruses": "Viruses"
    }
    return node['metadata']['name'].split(' ')[0] in kingdom_names

def process_taxonomy_data(data):
    if not isinstance(data, dict) or 'results' not in data:
        return None
    
    kingdom_counts = defaultdict(int)
    total_count = data.get('count', 0)
    
    # Process each result to find kingdom-level classifications
    for entry in data['results']:
        if 'metadata' in entry:
            name = entry['metadata']['name'].split(' ')[0]
            if name in ['Bacteria', 'Archaea', 'Eukaryota']:
                kingdom_counts[name] += 1
            # Special handling for viruses which might appear under different names
            elif 'vir' in name.lower():
                kingdom_counts['Viruses'] += 1
    
    # If we found no kingdoms but have results, try to count based on parent relationships
    if not kingdom_counts and data['results']:
        for entry in data['results']:
            if 'metadata' in entry and entry['metadata'].get('parent') == '131567':  # cellular organisms
                name = entry['metadata']['name'].split(' ')[0]
                if name in ['Bacteria', 'Archaea', 'Eukaryota']:
                    kingdom_counts[name] += 1
    
    if not kingdom_counts:
        return None
        
    # Calculate percentages
    formatted_results = {}
    for kingdom, count in kingdom_counts.items():
        percentage = (count / total_count) * 100
        formatted_results[kingdom] = {
            'count': count,
            'percentage': f"{percentage:.2f}%",
            'raw_count': f"{count:,}"
        }
    
    return formatted_results

def main():
    parser = argparse.ArgumentParser(description="Fetch Pfam taxonomic data")
    parser.add_argument("--input", required=True, help="File containing Pfam IDs, one per line")
    parser.add_argument("--output", help="Output JSON file (optional)")
    parser.add_argument("--debug", action="store_true", help="Print debug information")
    args = parser.parse_args()
    
    results = {}
    
    with open(args.input, 'r') as f:
        pfam_ids = [line.strip() for line in f]
    
    for pfam_id in pfam_ids:
        print(f"\nProcessing Pfam ID: {pfam_id}")
        data = fetch_pfam_data(pfam_id)
        
        if isinstance(data, str):
            print(f"  {data}")
            continue
        
        if args.debug:
            print(f"\nDEBUG: Found {data.get('count', 0)} total entries")
            if data['results']:
                print("\nDEBUG: Example taxonomic structure:")
                print(json.dumps(data['results'][0], indent=2))
        
        taxonomy = process_taxonomy_data(data)
        if taxonomy:
            results[pfam_id] = taxonomy
            print(f"Kingdom-level taxonomic breakdown (Total entries: {data.get('count', 0):,}):")
            for kingdom, stats in taxonomy.items():
                print(f"  {kingdom}:")
                print(f"    Count: {stats['raw_count']}")
                print(f"    Percentage: {stats['percentage']}")
        else:
            print("  No kingdom-level taxonomic data found")
        
        print("=" * 50)
    
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
            print(f"\nResults saved to {args.output}")

if __name__ == "__main__":
    main()