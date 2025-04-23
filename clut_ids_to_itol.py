import os
import csv
import argparse

def generate_itol_dataset(directory_path):
    """
    Generate an ITOL dataset file from text files in a directory.

    Parameters:
    directory_path (str): Path to the directory containing the text files.

    Returns:
    str: The generated ITOL dataset file content.
    """
    # Collect all unique IDs from all files
    all_ids = set()
    domain_data = {}

    # Gather all unique IDs and map domain files with assigned values
    for filename in os.listdir(directory_path):
        if filename.endswith(".txt"):
            file_path = os.path.join(directory_path, filename)
            domain = filename.replace('.txt', '')
            with open(file_path, "r") as file:
                ids_in_file = set(int(line.strip()) for line in file if line.strip())
                domain_data[domain] = ids_in_file
                all_ids.update(ids_in_file)

    # Sort IDs and domains for consistent output order
    all_ids = sorted(all_ids)
    domains = sorted(domain_data.keys())

    # Define colors based on values
    value_colors = {
        1: '#3F78C1',
        2: '#9163C1',
        3: '#44A043',
        4: '#D37538',
        5: '#B33D90',
        -1: '-1'
    }

    # Prepare dataset rows with presence (1-5) or absence (-1)
    dataset = []
    for id in all_ids:
        row_data = [id] + [1 if id in domain_data[domain] else -1 for domain in domains]
        dataset.append(row_data)

    # Assign colors to each field based on values
    field_colors = ['#000000']
    for domain in domains:
        color = value_colors.get(1)  # Default color for '1'
        field_colors.append(color)

    # Generate the ITOL dataset file content
    dataset_content = "\n".join([
        "DATASET_BINARY",
        "SEPARATOR\tCOMMA",
        "DATASET_LABEL,label1",
        "COLOR,#ff0000",
        "FIELD_SHAPES," + ",".join(['1'] * len(domains)),
        "FIELD_LABELS," + ",".join(domains),
        "FIELD_COLORS," + ",".join(field_colors),
        "DATA"] + [",".join(map(str, row)) for row in dataset])

    return dataset_content

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an ITOL dataset file from text files in a directory.")
    parser.add_argument("directory", help="Path to the directory containing the text files.")
    args = parser.parse_args()

    itol_dataset = generate_itol_dataset(args.directory)
    print(itol_dataset)
