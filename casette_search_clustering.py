import csv
import argparse
from collections import defaultdict
from tqdm import tqdm

def read_input_tsv(input_file):
    genes = {}
    with open(input_file, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in tqdm(reader, desc="Reading input file", unit=" genes"):
            if row.get('start') is None or row.get('end') is None or row.get('seqid') is None:
                continue
            gene_key = (row['seqid'], row['start'], row['end'])
            if gene_key not in genes or prioritize_source(row, genes[gene_key]):
                genes[gene_key] = row
    print(f"Read {len(genes)} genes from the input file after filtering incomplete entries.")
    return list(genes.values())

def prioritize_source(new_row, existing_row):
    return new_row.get('source') == 'img_core_v400' and existing_row.get('source') != 'img_core_v400'

def identify_clusters(genes, max_distance, min_pfams, core_pfams):
    clusters = []
    current_cluster = []
    cluster_id = 1
    sorted_genes = sorted(genes, key=lambda x: (x['seqid'], int(x['start'])))
    for i in tqdm(range(len(sorted_genes)), desc="Identifying clusters", unit=" genes"):
        if not current_cluster:
            current_cluster.append(sorted_genes[i])
            continue
        last_gene = current_cluster[-1]
        current_gene = sorted_genes[i]
        if last_gene['seqid'] == current_gene['seqid']:
            distance = int(current_gene['start']) - int(last_gene['end'])
            if distance <= max_distance:
                current_cluster.append(current_gene)
            else:
                if is_valid_cluster(current_cluster, min_pfams, core_pfams):
                    clusters.append((cluster_id, current_cluster))
                    cluster_id += 1
                current_cluster = [current_gene]
        else:
            if is_valid_cluster(current_cluster, min_pfams, core_pfams):
                clusters.append((cluster_id, current_cluster))
                cluster_id += 1
            current_cluster = [current_gene]
    if is_valid_cluster(current_cluster, min_pfams, core_pfams):
        clusters.append((cluster_id, current_cluster))
    print(f"Identified {len(clusters)} initial clusters.")
    return clusters

def is_valid_cluster(cluster, min_pfams, core_pfams):
    if len(cluster) < 3:
        return False
    unique_pfams = set()
    for gene in cluster:
        pfams = gene.get('pfamID', '').split(',')
        unique_pfams.update(pfams)
    if len(unique_pfams) < min_pfams:
        return False
    if not any(pfam in core_pfams for pfam in unique_pfams):
        return False
    return True

def calculate_cluster_score(cluster, core_pfams):
    unique_pfams = set()
    for gene in cluster:
        pfams = gene.get('pfamID', '').split(',')
        unique_pfams.update(pfams)
    num_core_pfams = len(unique_pfams.intersection(core_pfams))
    distances = [int(cluster[i+1]['start']) - int(cluster[i]['end']) for i in range(len(cluster)-1)]
    mean_distance = sum(distances) / len(distances) if distances else 0
    same_strand_percent = sum(1 for gene in cluster if gene.get('strand', '') == '+') / len(cluster)
    distance_factor = 1 / (mean_distance + 1) if mean_distance > 0 else 1
    score = num_core_pfams * 10 + distance_factor * 5 + same_strand_percent * 2
    return round(score, 1)

def merge_split_clusters(clusters, max_split_distance, core_pfams):
    merged_clusters = []
    i = 0
    while i < len(clusters):
        current_cluster = clusters[i]
        if i + 1 < len(clusters):
            next_cluster = clusters[i + 1]
            if (current_cluster[1][-1]['seqid'] == next_cluster[1][0]['seqid'] and
                int(next_cluster[1][0]['start']) - int(current_cluster[1][-1]['end']) <= max_split_distance and
                calculate_cluster_score(current_cluster[1], core_pfams) > 40 and
                calculate_cluster_score(next_cluster[1], core_pfams) > 40):
                merged_clusters.append((f"{current_cluster[0]}.1", current_cluster[1]))
                merged_clusters.append((f"{current_cluster[0]}.2", next_cluster[1]))
                i += 2
                continue
        merged_clusters.append(current_cluster)
        i += 1
    return merged_clusters

def write_output_files(clusters, output_file, summary_file, core_pfams):
    fieldnames = ['genomeID', 'geneID', 'pfamID', 'pfamDescription', 'seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    with open(output_file, 'w', newline='') as out_file, open(summary_file, 'w', newline='') as sum_file:
        out_writer = csv.writer(out_file, delimiter='\t')
        sum_writer = csv.writer(sum_file, delimiter='\t')
        out_writer.writerow(['clusterID', 'clusterScore'] + fieldnames)
        sum_writer.writerow(['clusterID', 'clusterScore', 'genomeID', 'firstGeneID', 'lastGeneID', 'clusterStart', 'clusterEnd', 'numUniquePfams', 'numCorePfams', 'meanDistance'])
        for cluster_id, cluster in tqdm(clusters, desc="Writing output files", unit=" clusters"):
            score = calculate_cluster_score(cluster, core_pfams)
            first_gene = cluster[0]
            last_gene = cluster[-1]
            unique_pfams = set()
            for gene in cluster:
                pfams = gene.get('pfamID', '').split(',')
                unique_pfams.update(pfams)
            num_unique_pfams = len(unique_pfams)
            num_core_pfams = len(unique_pfams.intersection(core_pfams))
            mean_distance = sum(int(cluster[i+1]['start']) - int(cluster[i]['end']) for i in range(len(cluster)-1)) / (len(cluster) - 1)
            for gene in cluster:
                pfam_domains = gene.get('pfamID', '').split(',')
                pfam_descriptions = gene.get('pfamDescription', '').split(',')
                for pfam_id, pfam_desc in zip(pfam_domains, pfam_descriptions):
                    out_writer.writerow([cluster_id, f"{score:.1f}", gene.get('genomeID', ''), gene.get('geneID', ''), 
                                         pfam_id, pfam_desc] + [gene.get(col, '') for col in fieldnames[6:]])
            sum_writer.writerow([cluster_id, f"{score:.1f}", first_gene.get('genomeID', ''), first_gene.get('geneID', ''), last_gene.get('geneID', ''),
                                 first_gene['start'], last_gene['end'], num_unique_pfams, num_core_pfams, f"{mean_distance:.1f}"])
    print(f"Wrote cluster details to {output_file}")
    print(f"Wrote cluster summaries to {summary_file}")

def main():
    parser = argparse.ArgumentParser(description='Identify gene clusters from TSV data.')
    parser.add_argument('--input-tsv', required=True, help='Input TSV file with gene data.')
    parser.add_argument('--max-distance', type=int, required=True, help='Maximum distance between neighboring genes.')
    parser.add_argument('--output-tsv', required=True, help='Output TSV file for clustered genes.')
    parser.add_argument('--summary-tsv', required=True, help='Output TSV file for cluster summaries.')
    parser.add_argument('--core-pfams', required=True, help='Comma-separated list of core Pfam IDs.')
    parser.add_argument('--min-pfams', type=int, required=True, help='Minimum number of unique Pfams in a cluster.')
    parser.add_argument('--max-split-distance', type=int, default=10000, help='Maximum distance for split clusters.')
    args = parser.parse_args()

    print("Starting gene cluster identification process...")
    print(f"Input file: {args.input_tsv}")
    print(f"Maximum distance between genes: {args.max_distance}")
    print(f"Minimum unique Pfams per cluster: {args.min_pfams}")
    print(f"Core Pfams: {args.core_pfams}")
    print(f"Maximum split cluster distance: {args.max_split_distance}")

    core_pfams = set(args.core_pfams.split(','))

    print("Reading input TSV file...")
    genes = read_input_tsv(args.input_tsv)

    print("Identifying clusters...")
    clusters = identify_clusters(genes, args.max_distance, args.min_pfams, core_pfams)

    print("Merging split clusters...")
    merged_clusters = merge_split_clusters(clusters, args.max_split_distance, core_pfams)

    print("Writing output files...")
    write_output_files(merged_clusters, args.output_tsv, args.summary_tsv, core_pfams)

    print("Process completed successfully.")

if __name__ == '__main__':
    main()