#!/usr/bin/env python

import argparse
import os
import re

def main(input_dir, output_dir):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Loop through all GFF files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".gff"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename)
            
            # Open the input and output files
            with open(input_path, "r") as input_file, open(output_path, "w") as output_file:
                # Loop through each line in the input file
                for line in input_file:
                    # Skip lines containing "gene"
                    if "\tgene\t" in line:
                        continue
                    # Replace the ID and locus_tag fields
                    line = re.sub(r'ID=product_(\w+);Parent=gene_(\w+);locus_tag=(\w+);', r'ID=cds_\3;locus_tag=\3;', line)
                    # Write the modified line to the output file
                    output_file.write(line)

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process GFF files")
    parser.add_argument("input_dir", help="path to input directory")
    parser.add_argument("output_dir", help="path to output directory")
    args = parser.parse_args()
    
    main(args.input_dir, args.output_dir)



print("Done")