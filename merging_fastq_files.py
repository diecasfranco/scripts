#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import gzip
import shutil

def merge_fastq(folderA, folderB, folderC):
    """
    Merges FASTQ or FASTQ.GZ files with the same name from two input folders into an output folder.
    """
    # Create output folder if it doesn't exist
    os.makedirs(folderC, exist_ok=True)
    
    # Loop through files in folderA
    for filename in os.listdir(folderA):
        if filename.endswith(".fastq") or filename.endswith(".fastq.gz"):
            fileA = os.path.join(folderA, filename)
            fileB = os.path.join(folderB, filename)
            output_file = os.path.join(folderC, filename)
            
            if os.path.exists(fileB):
                # Determine if the files are gzipped
                is_gzipped = filename.endswith(".gz")
                
                # Open the correct file types based on compression
                open_func = gzip.open if is_gzipped else open
                mode = "wb" if is_gzipped else "wb"  # Changed "w" to "wb" for non-gzipped files
                
                with open_func(output_file, mode) as out_f:
                    for f in [fileA, fileB]:
                        with open_func(f, "rb") as in_f:
                            shutil.copyfileobj(in_f, out_f)
                
                print(f"Merged: {filename}")
            else:
                print(f"Skipping: {filename} (no match in {folderB})")
    
    print("Merging complete!")

if __name__ == "__main__":  # Fixed the double asterisks
    parser = argparse.ArgumentParser(description="Merge FASTQ or FASTQ.GZ files with the same name from two folders.")
    parser.add_argument("-input_folder_1", required=True, help="Path to the first input folder")
    parser.add_argument("-input_folder_2", required=True, help="Path to the second input folder")
    parser.add_argument("-output_folder", required=True, help="Path to the output folder")
    
    args = parser.parse_args()
    merge_fastq(args.input_folder_1, args.input_folder_2, args.output_folder)
