#!/usr/bin/env python3

import re
import argparse
import os
import subprocess
import sys
from Bio import SeqIO

def read_patterns(pattern_file):
    patterns = []
    for record in SeqIO.parse(pattern_file, "fasta"):
        pattern = str(record.seq).replace('M', '.') # Replace 'M' with '.' (regex for any character)
        patterns.append(pattern)
    return patterns

def remove_pattern(sequence, patterns):
    for pattern in patterns:
        match = re.search(pattern, sequence)
        if match:
            return sequence[match.end():]
    return sequence

def process_fasta(fasta_file, patterns, db_file=None, show_progress=False):
    output_file = os.path.splitext(fasta_file)[0] + "_processed.fasta"
    with open(output_file, 'w') as output_handle:
        output_records = []
        total_records = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            total_records += 1
            original_seq = str(record.seq)
            new_seq = remove_pattern(original_seq, patterns)
            record.seq = new_seq
            output_records.append(record)
            if show_progress and total_records % 1000 == 0:
                print(f"Processed {total_records} records...", file=sys.stderr)

        SeqIO.write(output_records, output_handle, 'fasta')

    # Execute the usearch command on the output file if a database file is provided
    if db_file:
        usearch_cmd = f"usearch -orient {output_file} -db {db_file} -fastaout {os.path.splitext(output_file)[0]}_final.fasta -threads 20"
        subprocess.run(usearch_cmd, shell=True, check=True)

    return output_records

def process_directory(directory, patterns, db_file=None, show_progress=False, convert_fastq=False):
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            input_file = os.path.join(directory, filename)
            process_fasta(input_file, patterns, db_file, show_progress)
        elif filename.endswith(".fastq") and convert_fastq:
            sample_name = os.path.splitext(filename)[0]
            input_file = os.path.join(directory, filename)
            output_file = os.path.join(directory, f"{sample_name}.fasta")
            convert_to_fasta(input_file, output_file, sample_name)
            process_fasta(output_file, patterns, db_file, show_progress)

def convert_to_fasta(fastq_file, fasta_file, sample_name):
    cmd = f"vsearch -fastq_filter {fastq_file} -fastaout {fasta_file} -relabel {sample_name}. -fasta_width 0 --fastq_qmax 90"
    subprocess.run(cmd, shell=True, check=True)

def calculate_statistics(records):
    lengths = [len(record.seq) for record in records]
    min_length = min(lengths) if lengths else 0
    max_length = max(lengths) if lengths else 0
    avg_length = sum(lengths) / len(lengths) if lengths else 0
    return min_length, max_length, avg_length

def main():
    parser = argparse.ArgumentParser(description="Remove patterns from FASTA/FASTQ sequences and calculate statistics.")
    parser.add_argument('-i', '--input', required=True, help="Input FASTA or FASTQ file or directory containing FASTA/FASTQ files")
    parser.add_argument('-p', '--primers', required=True, help="FASTA file containing patterns")
    parser.add_argument('-db', '--db', help="Database file for usearch (optional)")
    parser.add_argument('-s', '--suffix', default="_trim", help="Suffix for output files when processing a directory")
    parser.add_argument('--stats', action='store_true', help="Show statistics for the FASTA/FASTQ file(s)")
    parser.add_argument('--progress', action='store_true', help="Show progress during processing")
    parser.add_argument('--convert-fastq', action='store_true', help="Convert FASTQ files to FASTA")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: {args.input} does not exist.")
        sys.exit(1)

    patterns = read_patterns(args.primers)
    if os.path.isdir(args.input):
        process_directory(args.input, patterns, args.db, args.progress, args.convert_fastq)
    else:
        if args.input.endswith(".fastq") and args.convert_fastq:
            sample_name = os.path.splitext(os.path.basename(args.input))[0]
            output_file = os.path.splitext(args.input)[0] + ".fasta"
            convert_to_fasta(args.input, output_file, sample_name)
            process_fasta(output_file, patterns, args.db, args.progress)
        else:
            process_fasta(args.input, patterns, args.db, args.progress)

    if args.stats:
        if os.path.isdir(args.input):
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
            records = list(SeqIO.parse(args.input, "fasta"))
            min_len, max_len, avg_len = calculate_statistics(records)
            print(f"Statistics for {args.input}:")
            print(f"Min length: {min_len}")
            print(f"Max length: {max_len}")
            print(f"Average length: {avg_len:.2f}")

if __name__ == "__main__":
    main()