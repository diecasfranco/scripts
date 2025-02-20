#!/usr/bin/env python3
import argparse
import os
import glob
import questionary

def substitute_sequence(sequence):
    substitution_dict = {'E': 'A', 'F': 'G', 'Q': 'C', 'P': 'T'}
    return ''.join(substitution_dict.get(char, char) for char in sequence)

def process_fasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write(line)
            else:
                corrected_sequence = substitute_sequence(line.strip())
                outfile.write(corrected_sequence + '\n')

def get_fasta_statistics(input_file):
    sequence_lengths = []
    current_sequence = []

    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                if current_sequence:
                    sequence_lengths.append(len(''.join(current_sequence)))
                    current_sequence = []
            else:
                current_sequence.append(line.strip())

        if current_sequence:
            sequence_lengths.append(len(''.join(current_sequence)))

    num_sequences = len(sequence_lengths)
    min_length = min(sequence_lengths) if sequence_lengths else 0
    max_length = max(sequence_lengths) if sequence_lengths else 0
    avg_length = sum(sequence_lengths) / num_sequences if num_sequences > 0 else 0

    return num_sequences, min_length, max_length, avg_length

def main():
    parser = argparse.ArgumentParser(description='Substitute specific letters in FASTA sequences.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file or directory containing FASTA files')
    parser.add_argument('-o', '--output', help='Output FASTA file')
    parser.add_argument('-s', '--statistics', action='store_true', help='Show statistics of the FASTA file(s)')

    args = parser.parse_args()

    if args.statistics:
        if os.path.isdir(args.input):
            fasta_files = glob.glob(os.path.join(args.input, '*.fasta'))
            for fasta_file in fasta_files:
                num_sequences, min_length, max_length, avg_length = get_fasta_statistics(fasta_file)
                print(f"Statistics for {fasta_file}:")
                print(f"  Number of sequences: {num_sequences}")
                print(f"  Min length: {min_length}")
                print(f"  Max length: {max_length}")
                print(f"  Average length: {avg_length:.2f}")
        else:
            num_sequences, min_length, max_length, avg_length = get_fasta_statistics(args.input)
            print(f"Statistics for {args.input}:")
            print(f"  Number of sequences: {num_sequences}")
            print(f"  Min length: {min_length}")
            print(f"  Max length: {max_length}")
            print(f"  Average length: {avg_length:.2f}")
    else:
        if os.path.isdir(args.input):
            fasta_files = glob.glob(os.path.join(args.input, '*.fasta'))
            
            overwrite = questionary.confirm("Do you want to overwrite the original files?").ask()
            
            if not overwrite:
                suffix = questionary.text("Please specify the suffix for new files (e.g., '_corrected'):").ask()
                for fasta_file in fasta_files:
                    output_file = f"{os.path.splitext(fasta_file)[0]}{suffix}.fasta"
                    process_fasta(fasta_file, output_file)
                    print(f"Substitutions have been made. Check the output file: {output_file}")
            else:
                for fasta_file in fasta_files:
                    process_fasta(fasta_file, fasta_file)
                    print(f"Substitutions have been made. The file has been overwritten: {fasta_file}")
        else:
            if args.output:
                process_fasta(args.input, args.output)
                print(f"Substitutions have been made. Check the output file: {args.output}")
            else:
                output_file = questionary.text("Please specify the output file name:").ask()
                process_fasta(args.input, output_file)
                print(f"Substitutions have been made. Check the output file: {output_file}")

if __name__ == '__main__':
    main()
