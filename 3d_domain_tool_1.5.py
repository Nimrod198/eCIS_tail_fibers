#!/usr/bin/env python3
import argparse
import logging
import numpy as np
from Bio.PDB import PDBParser, DSSP
from collections import defaultdict
import sys
import os
import glob

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Domain boundary analysis on a single structure from an AlphaFold PDB file using cross interactions only."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb_file", help="Path to the PDB file.")
    group.add_argument("--pdb_batch", help="Path to directory containing PDB files for batch processing.")
    parser.add_argument("--distance_threshold", type=float, default=8.0,
                        help="Distance threshold (in ?) for considering residue pairing (default: 8.0 ?).")
    parser.add_argument("--debug", action="store_true", help="Enable debug output for detailed reporting.")
    parser.add_argument("--min_domain_size", type=int, default=25,
                        help="Minimum number of residues for a domain to be considered significant (default: 25).")
    parser.add_argument("--max_subdomain_size", type=int, default=100,
                        help="Maximum number of residues for a subdomain (default: 100).")
    parser.add_argument("--overlap_threshold", type=float, default=0.8,
                        help="Overlap threshold (0-1) for merging similar subdomains (default: 0.8).")
    parser.add_argument("--output_tsv", type=str, help="Path to output TSV file (optional).")
    return parser.parse_args()

def init_logging(debug=False):
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(level=log_level,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        stream=sys.stdout)

def extract_secondary_structure(pdb_file, debug=False):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file)
    sec_structures = {}
    counts = {"H": 0, "B": 0, "L": 0}
    for key in dssp.keys():
        chain_id = key[0]
        resseq = key[1][1]
        dssp_code = dssp[key][2].strip()
        if dssp_code in ["H", "G", "I"]:
            ss = "H"
        elif dssp_code in ["E", "B"]:
            ss = "B"
        else:
            ss = "L"
        sec_structures[(chain_id, resseq)] = ss
        counts[ss] += 1
    if debug:
        logging.debug("Secondary structure summary from DSSP:")
        logging.debug(f"  Helices (H): {counts['H']}")
        logging.debug(f"  Sheets  (B): {counts['B']}")
        logging.debug(f"  Loops   (L): {counts['L']}")
    return sec_structures, model

def group_secondary_structures(model, sec_structures, debug=False):
    residue_to_unit = {}
    units = {}
    unit_id_counter = 1
    for chain in model:
        residues = [res for res in chain if res.has_id("CA") and (chain.id, res.id[1]) in sec_structures]
        residues.sort(key=lambda res: res.id[1])
        current_unit = []
        current_type = None
        for res in residues:
            key = (chain.id, res.id[1])
            ss = sec_structures.get(key, "L")
            if ss not in ["H", "B"]:
                if current_unit:
                    uid = f"unit_{unit_id_counter}"
                    unit_id_counter += 1
                    for r_key in current_unit:
                        residue_to_unit[r_key] = uid
                    units[uid] = {"chain": chain.id, "type": current_type,
                                  "residues": current_unit,
                                  "start": current_unit[0][1],
                                  "end": current_unit[-1][1]}
                    current_unit = []
                    current_type = None
                continue
            if current_unit and ss == current_type and (res.id[1] == current_unit[-1][1] + 1):
                current_unit.append(key)
            else:
                if current_unit:
                    uid = f"unit_{unit_id_counter}"
                    unit_id_counter += 1
                    for r_key in current_unit:
                        residue_to_unit[r_key] = uid
                    units[uid] = {"chain": chain.id, "type": current_type,
                                  "residues": current_unit,
                                  "start": current_unit[0][1],
                                  "end": current_unit[-1][1]}
                current_unit = [key]
                current_type = ss
        if current_unit:
            uid = f"unit_{unit_id_counter}"
            unit_id_counter += 1
            for r_key in current_unit:
                residue_to_unit[r_key] = uid
            units[uid] = {"chain": chain.id, "type": current_type,
                          "residues": current_unit,
                          "start": current_unit[0][1],
                          "end": current_unit[-1][1]}
    if debug:
        logging.debug("Secondary structure units formed:")
        for uid, details in units.items():
            logging.debug(f"  {uid}: Chain {details['chain']}, Type {details['type']}, Residues {details['start']}-{details['end']} (Count: {len(details['residues'])})")
    return residue_to_unit, units

def calculate_closest_residues(model, sec_structures, residue_to_unit, distance_threshold, debug=False):
    ca_coords = {}
    for chain in model:
        for residue in chain:
            if residue.has_id("CA"):
                key = (chain.id, residue.id[1])
                ca_coords[key] = residue["CA"].coord
    if debug:
        logging.debug(f"Total CA atoms considered: {len(ca_coords)}")

    closest_residues = {}
    pairing_details = {}

    for res_key, coord in ca_coords.items():
        if res_key not in residue_to_unit:
            continue
        min_dist = float("inf")
        closest_partner = None
        for other_key, other_coord in ca_coords.items():
            if res_key == other_key:
                continue
            if other_key not in residue_to_unit:
                continue
            if residue_to_unit[res_key] == residue_to_unit[other_key]:
                continue
            dist = np.linalg.norm(coord - other_coord)
            if dist < min_dist and dist <= distance_threshold:
                min_dist = dist
                closest_partner = other_key
        closest_residues[res_key] = closest_partner
        pairing_details[res_key] = (closest_partner, min_dist if closest_partner is not None else None)
        if debug and closest_partner:
            logging.debug(f"Residue {res_key} (Unit {residue_to_unit[res_key]}) paired with {closest_partner} (Unit {residue_to_unit[closest_partner]}) at {min_dist:.2f} ?")
    return closest_residues, pairing_details

def pair_secondary_structures(closest_residues, residue_to_unit, debug=False):
    pair_counts = defaultdict(int)
    for res_key, partner in closest_residues.items():
        if partner is not None:
            unit1 = residue_to_unit[res_key]
            unit2 = residue_to_unit[partner]
            if unit1 != unit2:
                pair = tuple(sorted([unit1, unit2]))
                pair_counts[pair] += 1
                if debug:
                    logging.debug(f"Recording cross interaction between {unit1} and {unit2} (current count: {pair_counts[pair]})")
    if debug:
        logging.debug(f"Total unique cross-unit interactions: {len(pair_counts)}")
    return pair_counts

def build_unit_graph(pair_counts, debug=False):
    graph = {}
    for (u1, u2) in pair_counts.keys():
        if u1 not in graph:
            graph[u1] = set()
        if u2 not in graph:
            graph[u2] = set()
        graph[u1].add(u2)
        graph[u2].add(u1)
    if debug:
        logging.debug("Unit graph connectivity:")
        for unit, neighbors in graph.items():
            logging.debug(f"  {unit}: {neighbors}")
    return graph

def connected_components(graph, debug=False):
    visited = set()
    components = []
    for node in graph:
        if node not in visited:
            stack = [node]
            comp = set()
            while stack:
                current = stack.pop()
                if current not in visited:
                    visited.add(current)
                    comp.add(current)
                    stack.extend(graph[current] - visited)
            components.append(comp)
    if debug:
        logging.debug("Connected components (domains) identified:")
        for comp in components:
            logging.debug(f"  Domain: {comp}")
    return components

def identify_subdomains(domain, units, max_subdomain_size):
    subdomains = []
    current_subdomain = set()
    current_size = 0
    for unit in sorted(domain, key=lambda u: units[u]["start"]):
        unit_size = len(units[unit]['residues'])
        if current_size + unit_size > max_subdomain_size and current_subdomain:
            subdomains.append(current_subdomain)
            current_subdomain = {unit}
            current_size = unit_size
        else:
            current_subdomain.add(unit)
            current_size += unit_size
    if current_subdomain:
        subdomains.append(current_subdomain)
    return subdomains

def get_domain_boundaries(domain_units, units):
    chain_ids = set()
    start_vals = []
    end_vals = []
    for uid in domain_units:
        unit = units.get(uid)
        if unit:
            chain_ids.add(unit["chain"])
            start_vals.append(unit["start"])
            end_vals.append(unit["end"])
    return chain_ids, min(start_vals), max(end_vals)

def merge_overlapping_subdomains(subdomains, units, overlap_threshold=0.9):
    filtered = []
    def get_boundaries(subdomain):
        _, start, end = get_domain_boundaries(subdomain, units)
        return start, end

    for sub in subdomains:
        merged = False
        s1, e1 = get_boundaries(sub)
        for i, fsub in enumerate(filtered):
            s2, e2 = get_boundaries(fsub)
            if s2 <= e1 and s1 <= e2:
                intersect = min(e1, e2) - max(s1, s2) + 1
                union = max(e1, e2) - min(s1, s2) + 1
                overlap = intersect / union if union > 0 else 0
                if overlap >= overlap_threshold:
                    filtered[i] = fsub.union(sub)
                    merged = True
                    break
        if not merged:
            filtered.append(sub)
    return filtered

def merge_and_identify_subdomains(domains, units, min_domain_size, max_subdomain_size, overlap_threshold):
    merged_domains = []
    for domain in domains:
        domain_size = sum(len(units[unit]['residues']) for unit in domain)
        if domain_size >= min_domain_size:
            subdomains = identify_subdomains(domain, units, max_subdomain_size)
            subdomains = [sub for sub in subdomains if sub != domain]
            filtered_subdomains = merge_overlapping_subdomains(subdomains, units, overlap_threshold)
            merged_domains.append((domain, filtered_subdomains))
    return merged_domains

def extract_b_scores(model, debug=False):
    b_scores = {}
    for chain in model:
        for residue in chain:
            if residue.has_id("CA"):
                key = (chain.id, residue.id[1])
                b_scores[key] = residue["CA"].get_bfactor()
    if debug:
        logging.debug(f"Extracted B-scores for {len(b_scores)} residues")
    return b_scores

def analyze_pdb_file(pdb_file, args, tsv_file=None):
    # Extract protein ID from the file name (before '__clust')
    pdb_basename = os.path.basename(pdb_file)
    protein_id = pdb_basename.split("__clust")[0] if "__clust" in pdb_basename else pdb_basename.split(".pdb")[0]
    
    logging.info(f"Analyzing file: {pdb_file}")
    sys.stdout.flush()

    try:
        sec_structures, model = extract_secondary_structure(pdb_file, debug=args.debug)
        total_ca = sum(1 for chain in model for res in chain if res.has_id("CA"))
        logging.info(f"Total residues with CA atoms in model: {total_ca}")

        residue_to_unit, units = group_secondary_structures(model, sec_structures, debug=args.debug)

        closest_residues, pairing_details = calculate_closest_residues(
            model, sec_structures, residue_to_unit, args.distance_threshold, debug=args.debug)
        paired_count = sum(1 for r, partner in closest_residues.items() if partner is not None)
        logging.info(f"Residues (in H/B units) with a valid cross interaction within {args.distance_threshold} ?: {paired_count}")

        pair_counts = pair_secondary_structures(closest_residues, residue_to_unit, debug=args.debug)
        unit_graph = build_unit_graph(pair_counts, debug=args.debug)
        domains = connected_components(unit_graph, debug=args.debug)

        merged_domains = merge_and_identify_subdomains(
            domains, units,
            args.min_domain_size,
            args.max_subdomain_size,
            args.overlap_threshold
        )

        # Extract B-scores (AlphaFold pLDDT scoring) from model CA ions.
        b_scores = extract_b_scores(model, debug=args.debug)
        
        domain_results = []
        for idx, (domain, subdomains) in enumerate(merged_domains, start=1):
            chain_ids, domain_start, domain_end = get_domain_boundaries(domain, units)
            # Compute average B-score for residues in the domain.
            domain_b_scores = []
            for uid in domain:
                for res_key in units[uid]["residues"]:
                    if res_key in b_scores:
                        domain_b_scores.append(b_scores[res_key])
            avg_b_score = np.mean(domain_b_scores) if domain_b_scores else float("nan")
            chain_str = ','.join(sorted(chain_ids))
            
            domain_results.append((protein_id, chain_str, domain_start, domain_end, avg_b_score))
            logging.info(f"Domain {idx} (Chain(s) {chain_str}): Residues {domain_start} - {domain_end}, Average B-score: {avg_b_score:.2f}")
        
        if tsv_file:
            for result in domain_results:
                protein_id, chain_str, domain_start, domain_end, avg_b_score = result
                tsv_file.write(f"{protein_id}\t{chain_str}\t{domain_start}\t{domain_end}\t{avg_b_score:.2f}\n")
        
        return domain_results
    except Exception as e:
        logging.error(f"Error processing {pdb_file}: {str(e)}")
        return []



def main():
    args = parse_arguments()
    init_logging(args.debug)

    if args.pdb_batch:
        # Batch processing mode
        directory = args.pdb_batch
        if not os.path.isdir(directory):
            logging.error(f"Directory not found: {directory}")
            sys.exit(1)
        
        # Create output TSV filename based on directory name
        dir_name = os.path.basename(os.path.normpath(directory))
        output_tsv = args.output_tsv if args.output_tsv else f"{dir_name}_dom_coordinates.tsv"
        
        # Find all PDB files in the directory
        pdb_files = glob.glob(os.path.join(directory, "*.pdb"))
        if not pdb_files:
            logging.error(f"No PDB files found in directory: {directory}")
            sys.exit(1)
        
        logging.info(f"Found {len(pdb_files)} PDB files in {directory}")
        
        # Process all PDB files and write results to a single TSV
        with open(output_tsv, "w") as tsv_out:
            tsv_out.write("ProteinID\tChains\tDomainStart\tDomainEnd\tAverageBScore\n")
            for pdb_file in pdb_files:
                analyze_pdb_file(pdb_file, args, tsv_out)
        
        logging.info(f"TSV output written to: {output_tsv}")
    
    else:
        # Single file processing mode
        pdb_file = args.pdb_file
        if not os.path.isfile(pdb_file):
            logging.error(f"File not found: {pdb_file}")
            sys.exit(1)
        
        # Extract protein ID from the file name
        pdb_basename = os.path.basename(pdb_file)
        protein_id = pdb_basename.split("__clust")[0] if "__clust" in pdb_basename else pdb_basename.split(".pdb")[0]
        
        # Create output TSV filename
        output_tsv = args.output_tsv if args.output_tsv else f"{protein_id}_dom_coordinates.tsv"
        
        with open(output_tsv, "w") as tsv_out:
            tsv_out.write("ProteinID\tChains\tDomainStart\tDomainEnd\tAverageBScore\n")
            analyze_pdb_file(pdb_file, args, tsv_out)
        
        logging.info(f"TSV output written to: {output_tsv}")

if __name__ == "__main__":
    main()
