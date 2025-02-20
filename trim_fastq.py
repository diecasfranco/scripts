#! /usr/bin/env python3

import os
import argparse
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from Bio import SeqIO

def process_file(file, input_folder, output_folder, min_length):
    if file.endswith(".fastq"):
        input_fastq = os.path.join(input_folder, file)
        output_fastq = os.path.join(output_folder, file)
        with open(output_fastq, 'w') as output_handle:
            for record in SeqIO.parse(input_fastq, "fastq"):
                if len(record.seq) >= min_length:
                    SeqIO.write(record, output_handle, "fastq")

def trim_fastq(input_folder, output_folder, min_length=250, num_cores=None):
    files = os.listdir(input_folder)
    num_cores = num_cores or cpu_count()
    with Pool(processes=num_cores) as pool:
        pool.starmap(process_file, [(file, input_folder, output_folder, min_length) for file in tqdm(files, desc="Processing files", unit="file")])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Trim FASTQ files based on minimum sequence length.')
    parser.add_argument('-i', '--input_folder', required=True, help='Path to the folder containing input FASTQ files.')
    parser.add_argument('-o', '--output_folder', required=True, help='Path to the folder where trimmed FASTQ files will be saved.')
    parser.add_argument('--min_length', type=int, default=250, help='Minimum sequence length to retain (default: 250)')
    parser.add_argument('--num_cores', type=int, default=None, help='Number of CPU cores to use for multiprocessing. Default is maximum available cores.')
    args = parser.parse_args()

    trim_fastq(args.input_folder, args.output_folder, args.min_length, args.num_cores)
