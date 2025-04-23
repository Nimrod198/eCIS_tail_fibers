import sys
import argparse
import random

# Enhanced Distinct Spring Palette
SPRING_PALETTE = [
    "#FF9A8B",  # Vibrant Coral
    "#4ECDC4",  # Teal
    "#45B7D1",  # Bright Sky Blue
    "#9B59B6",  # Deep Lilac
    "#F5D76E",  # Bright Sunflower Yellow
    "#E74C3C",  # Vivid Coral Red
    "#2ECC71",  # Emerald Green
    "#8E44AD"   # Rich Purple
]

def get_color(index):
    return SPRING_PALETTE[index % len(SPRING_PALETTE)]

def parse_input(input_text, color_map):
    result = []
    current_category = ""
    category_colors = {}
    categories = {}
    color_index = 0
    
    # First pass: identify all categories that share the first phrase
    for line in input_text.split('\n'):
        if line.startswith('>'):
            full_category = line.split('(')[0].strip()[1:]  # Remove '>' and trailing spaces
            parts = full_category.split('-')
            if parts[0] not in categories:
                categories[parts[0]] = []
            categories[parts[0]].append(full_category)

    # Second pass: process the data
    for line in input_text.split('\n'):
        if line.startswith('>'):
            full_category = line.split('(')[0].strip()[1:]  # Remove '>' and trailing spaces
            parts = full_category.split('-')
            if len(categories[parts[0]]) > 1:
                current_category = full_category
            else:
                current_category = parts[0]
            
            if current_category.lower() not in color_map and current_category not in category_colors:
                category_colors[current_category] = get_color(color_index)
                color_index += 1
        elif line and not line.startswith('>'):
            pfam = line.split()[1]
            color = color_map.get(current_category.lower(), category_colors.get(current_category))
            result.append(f"{current_category}\t{pfam}\t{color}")
    
    return result

def main():
    parser = argparse.ArgumentParser(description='Convert domain data to TSV format.')
    parser.add_argument('input_file', type=str, help='Input file path')
    parser.add_argument('output_file', type=str, help='Output file path')
    parser.add_argument('--colors', nargs='*', help='Color specifications (e.g., --colors viruses EFB036 bacteria FF0000)')

    args = parser.parse_args()

    color_map = {}
    if args.colors:
        for i in range(0, len(args.colors), 2):
            color_map[args.colors[i].lower()] = args.colors[i+1]

    with open(args.input_file, 'r') as f:
        input_text = f.read()

    result = parse_input(input_text, color_map)

    with open(args.output_file, 'w') as f:
        f.write("Category\tPfam\tColor\n")
        for line in result:
            f.write(line + '\n')

if __name__ == "__main__":
    main()
