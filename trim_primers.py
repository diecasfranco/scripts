#! /usr/bin/env python3

import re
import argparse
import os
from Bio import SeqIO

def read_patterns(pattern_file):
    patterns = []
    for record in SeqIO.parse(pattern_file, "fasta"):
        pattern = str(record.seq).replace('M', '.')  # Replace 'M' with '.' (regex for any character)
        patterns.append(pattern)
    return patterns

def remove_pattern(sequence, patterns):
    for pattern in patterns:
        match = re.search(pattern, sequence)
        if match:
            return sequence[match.end():]
    return sequence

def process_fasta(fasta_file, patterns):
    output_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        original_seq = str(record.seq)
        new_seq = remove_pattern(original_seq, patterns)
        record.seq = new_seq
        output_records.append(record)
    return output_records

def process_directory(directory, patterns, suffix):
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            input_file = os.path.join(directory, filename)
            output_file = os.path.join(directory, os.path.splitext(filename)[0] + suffix + ".fasta")
            processed_records = process_fasta(input_file, patterns)
            with open(output_file, 'w') as output_handle:
                SeqIO.write(processed_records, output_handle, 'fasta')

def calculate_statistics(records):
    lengths = [len(record.seq) for record in records]
    min_length = min(lengths) if lengths else 0
    max_length = max(lengths) if lengths else 0
    avg_length = sum(lengths) / len(lengths) if lengths else 0
    return min_length, max_length, avg_length

def main():
    parser = argparse.ArgumentParser(description="Remove patterns from FASTA sequences and calculate statistics.")
    parser.add_argument('-i', '--input', help="Input FASTA file or directory containing FASTA files")
    parser.add_argument('-o', '--output', help="Output FASTA file (ignored if input is a directory)")
    parser.add_argument('-p', '--primers', required=True, help="FASTA file containing patterns")
    parser.add_argument('-s', '--suffix', default="_trim", help="Suffix for output files when processing a directory")
    parser.add_argument('--stats', action='store_true', help="Show statistics for the FASTA file(s)")

    args = parser.parse_args()

    patterns = read_patterns(args.patterns)

    if os.path.isdir(args.input):
        process_directory(args.input, patterns, args.suffix)
        if args.stats:
            for filename in os.listdir(args.input):
                if filename.endswith(".fasta"):
                    input_file = os.path.join(args.input, filename)
                    records = list(SeqIO.parse(input_file, "fasta"))
                    min_len, max_len, avg_len = calculate_statistics(records)
                    print(f"Statistics for {filename}:")
                    print(f"Min length: {min_len}")
                    print(f"Max length: {max_len}")
                    print(f"Average length: {avg_len:.2f}\n")
    else:
        processed_records = process_fasta(args.input, patterns)
        with open(args.output, 'w') as output_handle:
            SeqIO.write(processed_records, output_handle, 'fasta')
        if args.stats:
            records = list(SeqIO.parse(args.input, "fasta"))
            min_len, max_len, avg_len = calculate_statistics(records)
            print(f"Statistics for {args.input}:")
            print(f"Min length: {min_len}")
            print(f"Max length: {max_len}")
            print(f"Average length: {avg_len:.2f}")

if __name__ == "__main__":
    main()
