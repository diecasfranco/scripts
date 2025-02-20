#! /usr/bin/env python3

import os
import pandas as pd

# get the input folder path from user
folderpath = input("Enter the path to the folder containing the fastq files: ")

# create an empty list to store the read counts
read_counts = []

# loop through each file in the folder
for filename in os.listdir(folderpath):
    # check if the file is a fastq file
    if filename.endswith('.fastq'):
        # open the fastq file
        with open(os.path.join(folderpath, filename), 'r') as file:
            # initialize a counter to keep track of the number of reads
            read_count = 0
            # loop through each line in the file
            for line in file:
                # if the line starts with '@', it's the header for a new read
                if line.startswith('@'):
                    read_count += 1
        # add the read count and filename to the list
        read_counts.append([filename, read_count])

# convert the list to a pandas dataframe
df = pd.DataFrame(read_counts, columns=['Filename', 'Read Count'])

# get the output table path from user
output_path = input("Enter the path and filename to save the output table (include .csv extension): ")

# save the dataframe to a csv file
df.to_csv(output_path, index=False)

print(f"Read counts for all fastq files in {folderpath} have been saved to {output_path}")
