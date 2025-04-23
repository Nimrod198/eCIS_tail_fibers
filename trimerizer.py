import os
import sys
import re

def process_fasta(input_file, output_file):
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            header = ''
            sequence = ''
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    if header and sequence:
                        outfile.write(f"{header}\n{sequence}:{sequence}:{sequence}\n")
                    header = line
                    sequence = ''
                else:
                    sequence += re.sub(r'[^A-Za-z]', '', line)
            
            if header and sequence:
                outfile.write(f"{header}\n{sequence}:{sequence}:{sequence}\n")
        print(f"Processed {input_file} -> {output_file}")
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")

def main(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    faa_files = [f for f in os.listdir(input_dir) if f.endswith('.faa')]
    print(f"Found {len(faa_files)} .faa files in {input_dir}")
    
    for filename in faa_files:
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, f"reformatted_{filename}")
        process_fasta(input_path, output_path)
    
    print(f"Processed files from {input_dir} and saved results to {output_dir}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_directory output_directory")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_dir = sys.argv[2]
    main(input_dir, output_dir)
