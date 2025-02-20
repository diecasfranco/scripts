#! /usr/bin/env python3
# Author: Piotr Łukasik, 1st Dec 2020

import sys

if len(sys.argv) != 3:
	sys.exit('This script adds up zOTU counts and adds the totals to the fasta file heading before OTU picking.\n'
           'This is necessary for the usearch OTU picking algorithm to function properly\n'
	         'Usage: ./add_values_to_zOTU_fasta.py <zOTU_count_table> <zOTU_fasta>\n')

Script, zOTU_counts, zOTU_fasta = sys.argv


CT_Table = open(zOTU_counts, "r")
zOTU_dict = {}

for line in CT_Table:
    if not line.startswith("#OTU"):
        LINE = line.strip().split()
        heading = LINE[0]
        total = 0
        for pos in range(1,len(LINE)):
            total += int(LINE[pos])
        zOTU_dict[heading] = total

CT_Table.close()

OLD_FASTA = open(zOTU_fasta, "r")
NEW_FASTA = open("new_" + zOTU_fasta, "w")

for line in OLD_FASTA:
    if line.startswith(">"):
        heading = line.strip().strip(">")
        print(">%s;size=%d" % (heading, zOTU_dict[heading]), file=NEW_FASTA)
    else:
        print(line.strip(), file=NEW_FASTA)

OLD_FASTA.close()
NEW_FASTA.close()
