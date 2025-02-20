#! /usr/bin/env python3

import argparse
import os

def extract_sequences(input_file, output_dir):
    sequences_by_identifier = {}

    with open(input_file, "r") as infile:
        for line in infile:
            line = line.rstrip()  # Remove newline character at the end of the line
            if line.startswith(">"):
                _, header_info = line.split(" ", 1)  # Split only once at the first space
                identifier = header_info.strip().replace(" ", "_")
                sequences_by_identifier.setdefault(identifier, []).append(line)
            elif identifier:
                sequences_by_identifier[identifier].append(line)
    for identifier, sequences in sequences_by_identifier.items():
        identifier = identifier.replace(" ", "_")  # Replace spaces with underscores
        identifier = identifier.replace("/", "-")  # Replace slashes with hyphens
        output_file = os.path.join(output_dir, f"{identifier}.faa")
        with open(output_file, "w") as outfile:
            outfile.write("\n".join(sequences))

def main():
    parser = argparse.ArgumentParser(description="Extract sequences with the same identifier from a FASTA file and save them to separate files.")
    parser.add_argument("-i", "--input_file", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output_dir", default="proteins_folder", help="Output directory for sequences with the same identifier")

    args = parser.parse_args()

    # Create output directory if it does not exist
    os.makedirs(args.output_dir, exist_ok=True)

    extract_sequences(args.input_file, args.output_dir)

if __name__ == "__main__":
    main()
