#! /usr/bin/env python3
import argparse

# Create command-line argument parser
parser = argparse.ArgumentParser(description="Update FASTA headers with taxonomy information")
parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
parser.add_argument("-t", "--taxonomy", required=True, help="Input taxonomy file")
args = parser.parse_args()

# Read the taxonomy file and create a dictionary to map codes to taxonomic labels
taxonomy_map = {}
with open(args.taxonomy, 'r') as taxonomy_file:
    for line in taxonomy_file:
        code, taxonomy = line.strip().split('\t')
        taxonomy_map[code] = taxonomy

# Read the FASTA file and create a new FASTA file with updated headers
with open(args.fasta, 'r') as input_fasta, open('output.fasta', 'w') as output_fasta:
    for line in input_fasta:
        if line.startswith('>'):
            code = line[1:].strip()  # Remove '>' and any leading/trailing whitespace
            if code in taxonomy_map:
                # Extract and format the taxonomy labels
                taxonomy = taxonomy_map[code]
                labels = taxonomy.split('; ')
                formatted_labels = ';'.join([f"{label.split('__')[-1]}" for label in labels])

                new_header = f'>{code};tax={formatted_labels}'
                output_fasta.write(new_header + '\n')
            else:
                # If the code is not found in the taxonomy file, keep the original header
                output_fasta.write(line)
        else:
            # Keep the sequence lines unchanged
            output_fasta.write(line)

print("FASTA headers have been updated and saved to output.fasta.")

