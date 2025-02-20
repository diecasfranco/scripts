#! /usr/bin/env python3

import os
from Bio import SeqIO
import argparse

def fasta_to_fastq(fasta_file, fastq_file, quality_score=40):
    """
    Convert a FASTA file to a FASTQ file with a default quality score.

    :param fasta_file: Input FASTA file path
    :param fastq_file: Output FASTQ file path
    :param quality_score: Quality score to assign to each base (default: 40)
    """
    with open(fasta_file, "r") as fasta_handle, open(fastq_file, "w") as fastq_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            # Assign a default quality score to each base
            record.letter_annotations["phred_quality"] = [quality_score] * len(record.seq)
            SeqIO.write(record, fastq_handle, "fastq")

def convert_all_fasta_in_folder(folder, quality_score=40):
    """
    Convert all FASTA files in the specified folder to FASTQ files.

    :param folder: Folder containing the FASTA files
    :param quality_score: Quality score to assign to each base (default: 40)
    """
    for filename in os.listdir(folder):
        if filename.endswith(".fasta"):
            fasta_file = os.path.join(folder, filename)
            fastq_file = os.path.join(folder, os.path.splitext(filename)[0] + ".fastq")
            fasta_to_fastq(fasta_file, fastq_file, quality_score)
            print(f"Converted {fasta_file} to {fastq_file} with quality score {quality_score}")

def main():
    parser = argparse.ArgumentParser(description="Convert all FASTA files in a folder to FASTQ files with default quality scores.")
    parser.add_argument("folder", help="Folder containing the FASTA files")
    parser.add_argument("--quality_score", type=int, default=40, help="Quality score to assign to each base (default: 40)")

    args = parser.parse_args()
    convert_all_fasta_in_folder(args.folder, args.quality_score)

if __name__ == "__main__":
    main()
