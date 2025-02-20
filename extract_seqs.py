#! /usr/bin/env python3

import argparse

def extract_sequences(input_file, output_file, protein_name):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        found_sequence = False
        for line in infile:
            if protein_name in line:
                found_sequence = True
                outfile.write(line)  # Write the header
            elif found_sequence and line.startswith(">"):
                found_sequence = False
            elif found_sequence:
                outfile.write(line.strip() + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract sequences matching a specified protein name from a FASTA file.")
    parser.add_argument("-i", "--input_file", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output_file", required=True, help="Output file for sequences matching the protein name")
    parser.add_argument("-name", "--protein_name", required=True, help="Protein name to search for in the FASTA file")

    args = parser.parse_args()

    extract_sequences(args.input_file, args.output_file, args.protein_name)

if __name__ == "__main__":
    main()


print("Done - Salud")