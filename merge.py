#! /usr/bin/env python3

import os, argparse, itertools, sys, multiprocessing, csv, time

# define command-line options
parser = argparse.ArgumentParser(description='This script is intended to produce a table with both contig information')
parser.add_argument("-c", "--coverage_table", metavar='<coverage_table>', help="the path to the coverage table file", required=True)
parser.add_argument("-b", "--blastn_table", metavar='<blastn_table>', help="the path to the blastn table file", required=True)
parser.add_argument("-o", "--output_path", metavar='<output_path>', help="the path to the output file", required=True)

args = parser.parse_args()
coverage_table = args.coverage_table
blastn_table = args.blastn_table
output_path = args.output_path

# import taxonomy table as a dictionary
with open("{}".format(blastn_table), "r") as Table:
    Table_tax = {}
    for line in Table: 
        LINE = line.strip().split()
        header = LINE[0]
        taxonomy = LINE[1]
        Table_tax[header] = taxonomy
with open("{}".format(coverage_table), "r") as Table:
    Table_cov = []
    for line in Table: 
        LINE = line.strip().split()
        Table_cov.append(LINE)

row_count = len(Table_cov)

# add the taxonomy from blastn table
Table_cov[0].append("blastn")
for row in range(1, row_count):
    if Table_cov[row][0] in Table_tax.keys(): 
        header = Table_cov[row][0]
        taxonomy = Table_tax[header]
        Table_cov[row].append(taxonomy)
    else: 
        Table_cov[row].append("NA")
col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "blastn"]
col_pos = []
for name in col_keep: 
    col_pos.append(Table_cov[0].index(name))

with open(output_path, "w") as OUTPUT:
    for row in Table_cov: 
        col_position = 0
        for col in row[:-1]:
            if col_position in col_pos:
                print(col, file = OUTPUT, end = "\t")
            col_position += 1
        print(row[-1], file = OUTPUT)

print("Done")

