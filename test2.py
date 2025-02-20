#! /usr/bin/env python3

import sys, os, re, subprocess


##### Setting names of output files
Output_table = "zotu_table_expanded.txt"

##### Setting up the key arrays --- LIST for keeping sequences in order, and DICT for managing sequence info
zOTU_list = []
zOTU_dict = {}


##### Opening zOTU table

COUNTS = open("all_libraries_zotu_table.txt", "r")

for line in COUNTS:
    if line.startswith("#"):
        COUNT_headings = line.strip().split()[1:]    ### Keeping the names of libraries
    else:
        LINE = line.strip().split()
        zOTU_list.append(LINE[0])
        zOTU_dict[LINE[0]] = [LINE[1:]]

COUNTS.close()


##### Adding taxonomy info to DICT

TAX = open("zotus.tax", "r")

for line in TAX:
    LINE = line.strip().split()
    if LINE[0] in zOTU_list:
        if len(LINE) > 1:
            zOTU_dict[LINE[0]].append(LINE[1])
        else:
            zOTU_dict[LINE[0]].append("unassigned")
    else:
        print('FATAL ERROR! Taxonomy file contains zOTUs that are not in zOTU count table! ---', LINE[0])
        sys.exit()

TAX.close()

##### Adding sequences from the FASTA file to DICT
FASTA = open("new_zotus.fasta", "r")
Sequence = ''
Seq_heading = FASTA.readline().strip().strip(">")

for line in FASTA:   # Copying the sequence (potentially spread across multiple lines) to a single line
    if line.startswith('>'):    # if the line contains the heading
        if Seq_heading not in zOTU_list and Seq_heading != "":     # EXIT if the previous Seq_heading is not in a list!
            print('FATAL ERROR! Fasta file contains zOTUs that are not in zOTU count table! ---', Seq_heading)
            sys.exit()
            
        zOTU_dict[Seq_heading].append(Sequence) # save the existing Seq_heading and Sequence to a DICT
        Sequence = ''    # clear sequence
        Seq_heading = line.strip().strip(">")  # use the current line as the new heading!

    else:
        Sequence = Sequence + line.strip().upper()

zOTU_dict[Seq_heading].append(Sequence) # Saves the final sequence (Seq_heading and Sequence) to a list

FASTA.close()

##### Adding zOTU - OTU relationship info to DICT

RELS = open("zotu_otu_relationships.txt", "r")

for line in RELS:
    LINE = line.strip().split()
    
    zOTU = re.search(r"^Zotu\d+", LINE[0])[0]
    if zOTU not in zOTU_list:
        print('FATAL ERROR! Relationship file contains zOTUs that are not in zOTU count table! --- ', zOTU)
        sys.exit()
    
    if LINE[1].startswith("otu"):
        zOTU_dict[zOTU].append(LINE[1])
    
    elif  LINE[1] == "noisy_chimera" or LINE[1] == "perfect_chimera" or LINE[1] == "match_chimera" or re.search("Chimera", LINE[2]) != None:
        zOTU_dict[zOTU].append("Chimera")

    elif (LINE[1] == "match" or LINE[1] == "perfect") and re.search("OTU\\d+", LINE[2]) != None:
        OTU_ID = re.search("OTU\\d+", LINE[2])[0].lower()
        zOTU_dict[zOTU].append(OTU_ID)


        
    else:
        print("Relationship file contains a term that I have not considered")
        sys.exit()

RELS.close()


##### Outputting the Expanded Count Table
OUTPUT_TABLE = open(Output_table, "w")

print("OTU_ID", "OTU_assignment", "Taxonomy", "Sequence", "Total", sep = "\t", end = "\t", file = OUTPUT_TABLE)
for item in COUNT_headings:
    print(item, end = "\t", file = OUTPUT_TABLE)
print("", file = OUTPUT_TABLE)

for zOTU in zOTU_list:
    Total = 0
    for no in zOTU_dict[zOTU][0]:
        Total += int(no)
    
    # Terms in DICT: 'Zotu32': [['0', '1', '100'], 'd:Bacteria(1.00)...', 'TACGT...', 'otu8']
    # I want to export: "OTU_ID", "OTU_assignment"[3], "Taxonomy"[1], "Sequence"[2], "Total"
    print(zOTU, zOTU_dict[zOTU][3], zOTU_dict[zOTU][1], zOTU_dict[zOTU][2], Total, sep = "\t", end = "\t", file = OUTPUT_TABLE)
    
    for no in zOTU_dict[zOTU][0]:
        print(no, end = "\t", file = OUTPUT_TABLE)
    
    print("", file = OUTPUT_TABLE)

OUTPUT_TABLE.close()

print("zOTU_Table_expanded ready")


### Creating OTU_Table:

OTU = open("zotu_table_expanded.txt", "r")
OTU_TABLE = []
for line in OTU:
    LINE = line.strip().split()
    if line.startswith("OTU_ID"):
        COUNT_headings = line.strip().split()[4:]    ### Keeping the names of libraries
    else:
        OTU_TABLE.append(LINE)   
OTU.close()

otu_dict = {}
for row_no in range(0, len(OTU_TABLE)):
    otu_key = OTU_TABLE[row_no][1]
    if not otu_key in otu_dict.keys():
         otu_dict[otu_key] = OTU_TABLE[row_no][4:]
    else:
        otu_dict[otu_key] = [sum(map(int, i)) for i in list(zip(otu_dict[otu_key], OTU_TABLE[row_no][4:]))]

##### Adding taxonomy info to DICT
TAX = open("otus.tax", "r")
OTU_TAX = []
for line in TAX:
    LINE = line.strip().split()
    OTU_TAX.append(LINE)

### Lowering the #OTU in Taxonomy file:
for list in OTU_TAX:
    list[0] = list[0].lower()       
        

for row_no in range(0, len(OTU_TAX)):   
    if OTU_TAX[row_no][0] in otu_dict.keys():
        if len(OTU_TAX[row_no]) > 1:
            otu_dict[OTU_TAX[row_no][0]].append(OTU_TAX[row_no][1])
        else:
            otu_dict[OTU_TAX[row_no][0]].append("unassigned")
TAX.close()


                
###We are adding 1 to the end of our dictionary to 
for row_no in range(0, len(OTU_TABLE)):
    if otu_dict[OTU_TABLE[row_no][1]][-1] != 1 and OTU_TABLE[row_no][1] in otu_dict.keys():
        otu_dict[OTU_TABLE[row_no][1]].append(OTU_TABLE[row_no][3])
        otu_dict[OTU_TABLE[row_no][1]].append(1)      

                
COUNT_headings.insert(0,"#OTU")
COUNT_headings.insert(1,"Taxonomy")
COUNT_headings.insert(2,"Sequence")
data = []
data.append(COUNT_headings)      
for otu in otu_dict.keys():
    data.append([otu] + [otu_dict[otu][-3]] + [otu_dict[otu][-2]] + otu_dict[otu][:-3])
 


import subprocess

# Write OTU table
with open("OTU_Table.txt", "w") as bigFile:
    for LINE in data:
        print("\t".join(map(str, LINE)), file=bigFile)

# Run sed command to remove text inside parentheses
subprocess.run(["sed", "-i", "-E", "s/\\([^()]*\\)//g", "OTU_Table.txt"])
subprocess.run(["sed", "-i", "-E", "s/\\([^()]*\\)//g", "zotu_table_expanded.txt"])
print("OTU_Table ready!")

ascii_logo = r"""                                             
 _____ _            __  __ _ _         ____                       
|_   _| |__   ___  |  \/  (_) |_ ___  |  _ \ ___   ___  _ __ ___  
  | | | '_ \ / _ \ | |\/| | | __/ _ \ | |_) / _ \ / _ \| '_ ` _ \ 
  | | | | | |  __/ | |  | | | ||  __/ |  _ < (_) | (_) | | | | | |
  |_| |_| |_|\___| |_|  |_|_|\__\___| |_| \_\___/ \___/|_| |_| |_|                                                               
😎 
                                                                 
"""
print(ascii_logo)
print("The Mite Room® Sante - goed bezig! Salud! Na zdrowie! Gānbēi (干杯)!")
