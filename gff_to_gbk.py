#!/usr/bin/env python

from Bio import SeqIO

def convert_gff_to_genbank(input_file, output_file):
    # Load the input file as a SeqRecord object
    record = SeqIO.read(input_file, "fasta")

    # Add the "molecule_type" annotation if missing
    if "molecule_type" not in record.annotations:
        record.annotations["molecule_type"] = "DNA"

    # Convert the input file to a GenBank file
    SeqIO.write(record, output_file, "genbank")

    # Print a success message
    print(f"Conversion complete. GenBank file saved as {output_file}.")

if __name__ == "__main__":
    import argparse

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Convert a GFF and associated FASTA file into GenBank format.")
    parser.add_argument("input_file", help="Path to the input file (GFF and FASTA combined)")
    parser.add_argument("output_file", help="Path to the output GenBank file")
    args = parser.parse_args()

    # Call the conversion function
    convert_gff_to_genbank(args.input_file, args.output_file)
