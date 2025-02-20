#!/usr/bin/env python3

import argparse
from Bio import SeqIO

# Create the parser
parser = argparse.ArgumentParser(description="Remove sequences from a FASTA file.")

# Add the arguments
parser.add_argument('input_fasta', type=str, help='The input FASTA file.')
parser.add_argument('output_fasta', type=str, help='The output FASTA file.')
parser.add_argument('remove_list', type=str, help='The list of sequence identifiers to remove.')

# Parse the arguments
args = parser.parse_args()

# Load the list of sequence identifiers to remove
with open(args.remove_list, 'r') as remove_file:
    remove_ids = set(line.strip() for line in remove_file)

# Open the input FASTA file and the output file
with open(args.input_fasta, 'r') as input_fasta, open(args.output_fasta, 'w') as output_fasta:
    for record in SeqIO.parse(input_fasta, 'fasta'):
        if record.id not in remove_ids:
            SeqIO.write(record, output_fasta, 'fasta')