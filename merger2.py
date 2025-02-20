#!/usr/bin/env python3

import os
import argparse
import datetime  # Added for getting the current date

# Define command-line options
parser = argparse.ArgumentParser(description='This script processes coverage and blastn tables to produce an output table with contig information.')
parser.add_argument("-c", "--coverage_table", metavar='<coverage_table>', help="Path to the coverage table file", required=True)
parser.add_argument("-b", "--blastn_table", metavar='<blastn_table>', help="Path to the blastn table file", required=True)
parser.add_argument("-o", "--output_path", metavar='<output_path>', help="Path to the output file", required=True)
parser.add_argument("--full_db", action="store_true", help="Use this flag if the blast table has the 16th and 17th columns in the format")
parser.add_argument("--full_name", action="store_true", help="Use this flag to print the whole scientific name provided by the blastn table")

args = parser.parse_args()
coverage_table = args.coverage_table
blastn_table = args.blastn_table
output_path = args.output_path
use_full_db_format = args.full_db
include_full_name = args.full_name

with open("{}".format(blastn_table), "r") as Table:
    Table_tax = {}
    for line in Table:
        LINE = line.strip().split("\t")  # Use tab as the delimiter
        header = LINE[0]
        if use_full_db_format:
            taxonomy = LINE[14] # Concatenate the 15th and 16th columns for taxonomy
            kingdom = LINE[15]  # Use the 16th column for kingdom
            full_name = LINE[16] if include_full_name else "NA"  # Use the last column for full name if --full_name is provided
        else:
            taxonomy = LINE[1]  # Use the 2nd column for taxonomy
            kingdom = "NA"  # Set to "NA" as it's not available in this format
            full_name = "NA"  # Set to "NA" as it's not available in this format
        Table_tax[header] = (taxonomy, kingdom, full_name)


with open("{}".format(coverage_table), "r") as Table:
    Table_cov = []
    for line in Table:
        LINE = line.strip().split()
        Table_cov.append(LINE)

row_count = len(Table_cov)

# Add the taxonomy, kingdom and full_name from blastn table
Table_cov[0].append("taxonomy")  # Header for the new column
Table_cov[0].append("kingdom")   # Header for the new column
if include_full_name:
    Table_cov[0].append("full_name")  # Header for the new column
for row in range(1, row_count):
    if Table_cov[row][0] in Table_tax.keys():
        header = Table_cov[row][0]
        taxonomy, kingdom, full_name = Table_tax[header]
        Table_cov[row].extend([taxonomy, kingdom, full_name if include_full_name else "NA"])
    else:
        Table_cov[row].extend(["NA", "NA", "NA" if include_full_name else ""])

col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "taxonomy", "kingdom"]
if include_full_name:
    col_keep.append("full_name")
col_pos = []
for name in col_keep:
    col_pos.append(Table_cov[0].index(name))

# Add the 15th and 16th columns if using the full database format
if use_full_db_format:
    col_pos.append(15)  # 16th column
    col_pos.append(16)  # 17th column

with open(output_path, "w") as OUTPUT:
    for row in Table_cov:
        col_position = 0
        for col in row[:-1]:
            if col_position in col_pos:
                print(col, file=OUTPUT, end="\t")
            col_position += 1
        print(row[-1], file=OUTPUT)

# Print the current date and a blinking eye emoticon at the end
current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print("\nScript completed on:", current_date)
print("Table is done! \U0001F609 \U0001F60E \U0001F37A")

print("Done - Salud")