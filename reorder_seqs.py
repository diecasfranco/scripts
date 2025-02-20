#! /usr/bin/env python3
### Script to reorder the sequences
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())

def calculate_similarity(seq1, seq2, length=50):
    """Calculate the similarity between two sequences over a given length."""
    return sum(a == b for a, b in zip(seq1[:length], seq2[:length])) / length

def is_reverse_complement(seq, reference_seq, threshold=0.8):
    """Determine if a sequence is a reverse complement of the reference sequence."""
    forward_similarity = calculate_similarity(seq, reference_seq)
    reverse_similarity = calculate_similarity(reverse_complement(seq), reference_seq)
    
    if reverse_similarity > forward_similarity and reverse_similarity > threshold:
        return True
    return False

def process_fasta_file(input_file, output_file, threshold):
    """Process a single FASTA file and write corrected sequences to the output file."""
    records = list(SeqIO.parse(input_file, "fasta"))

    if not records:
        print(f"No sequences found in the input file {input_file}.")
        return

    reference_seq = str(records[0].seq)
    corrected_records = [records[0]]  # Keep the first record as is

    for record in records[1:]:
        seq = str(record.seq)
        if is_reverse_complement(seq, reference_seq, threshold):
            record.seq = Seq(reverse_complement(seq))
            record.description += " [Reverse Complemented]"
        corrected_records.append(record)

    # Write the corrected sequences to a new FASTA file
    SeqIO.write(corrected_records, output_file, "fasta")
    print(f"Corrected sequences have been written to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Process all FASTA files in a folder, identify reverse complements, and reorder sequences.",
        epilog="Credits: DCF, Claude and chatGPT"
    )
    parser.add_argument("input_folder", help="Path to the input folder containing FASTA files")
    parser.add_argument("--output_folder", help="Path to the output folder for corrected FASTA files", default=None)
    parser.add_argument("--threshold", type=float, default=0.8, help="Threshold for reverse complement similarity (default: 0.8)")

    args = parser.parse_args()

    input_folder = args.input_folder
    output_folder = args.output_folder if args.output_folder else input_folder
    threshold = args.threshold

    if not os.path.isdir(input_folder):
        print(f"The specified input folder does not exist: {input_folder}")
        return

    if not os.path.isdir(output_folder):
        print(f"The specified output folder does not exist: {output_folder}")
        return

    fasta_files = [f for f in os.listdir(input_folder) if f.endswith(".fasta")]

    if not fasta_files:
        print(f"No FASTA files found in the specified folder: {input_folder}")
        return

    for fasta_file in fasta_files:
        input_file = os.path.join(input_folder, fasta_file)
        base_name = os.path.splitext(fasta_file)[0]
        output_file = os.path.join(output_folder, f"{base_name}_reorder.fasta")
        process_fasta_file(input_file, output_file, threshold)

if __name__ == "__main__":
    main()
