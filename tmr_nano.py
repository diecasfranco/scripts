#!/usr/bin/env python3

import sys
import os
import re
import argparse
import subprocess
from tqdm import tqdm
from datetime import datetime

# Argument parser
parser = argparse.ArgumentParser(description="Cluster OTUs and assign taxonomy using USEARCH and VSEARCH.")
parser.add_argument("-db", "--db_path", help="Path to the database file")
parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads (default: 8)")

args = parser.parse_args()
db_path = args.db_path
j = args.threads

# Get current time
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(f"Current date and time: {current_time}")

# Copy and rename consensus sequence files
os.system(r"""find . -name "*_consensussequences.fasta" -exec cp "{}" . \;""")
for filename in os.listdir():
    if filename.endswith(".fasta"):
        new_name = re.sub(r"_consensussequences", "", filename)
        if new_name != filename:
            os.rename(filename, new_name)

# Modify fasta headers
counter = 0
os.system('clear')

for file_name in os.listdir('.'):
    if file_name.endswith('.fasta'):
        SampleName = file_name.replace('.fasta', '')
        with open(file_name, 'r') as input_file:
            lines = input_file.readlines()

        with open(file_name, 'w') as output_file:
            for line in lines:
                if line.startswith('>'):
                    line = re.sub(rf'>consensus_{SampleName}_[0-9_]+\((\d+)\)', 
                                  rf'>{SampleName}._{counter};size=\1', line)
                    counter += 1
                output_file.write(line)

# Dereplication
os.system("""for file in *.fasta; do SampleName=`basename $file .fasta`; vsearch -sortbysize "$SampleName".fasta --output "$SampleName".sorted.fasta -minsize 2; done""")

# Denoising
print("Denoising..................... ", end="")
os.system("""for file in *sorted.fasta; do SampleName=`basename $file .sorted.fasta`; usearch -unoise3 "$SampleName".sorted.fasta -zotus "$SampleName".zotus.fasta -tabbedout "$SampleName".denoising.summary.txt -minsize 1; done""")
os.system("mkdir sorted && mv *sorted.fasta sorted/")
os.system("mkdir denoising_summary && mv *denoising.summary.txt denoising_summary/")

# OTU Table
os.system(f"""for file in *zotus.fasta; do SampleName=`basename $file .zotus.fasta`; usearch -otutab "$SampleName".fasta -zotus "$SampleName".zotus.fasta -otutabout "$SampleName"_zotu_table.txt -threads 20; done""")

# === Native Python version of add_seq_to_zotu ===

def add_seq_to_zotu(zotu_counts, zotu_fasta, output):
    zOTU_list = []
    zOTU_dict = {}

    with open(zotu_counts, 'r') as counts:
        for line in counts:
            if line.startswith("#"):
                COUNT_headings = line.strip().split()[2:]
            else:
                parts = line.strip().split()
                zOTU_list.append(parts[0])
                zOTU_dict[parts[0]] = [parts[1:]]

    with open(zotu_fasta, 'r') as fasta:
        sequence = ''
        heading = fasta.readline().strip().strip(">")

        for line in fasta:
            if line.startswith(">"):
                if heading not in zOTU_list and heading != "":
                    print(f"FATAL ERROR! Fasta file contains zOTUs not in count table: {heading}")
                    sys.exit(1)
                zOTU_dict[heading].append(sequence)
                sequence = ''
                heading = line.strip().strip(">")
            else:
                sequence += line.strip().upper()

        # Final sequence
        zOTU_dict[heading].append(sequence)

    with open(output, 'w') as out:
        out.write("OTU_ID\tSequence\t" + "\t".join(COUNT_headings) + "\n")
        for zotu in zOTU_list:
            out.write(f"{zotu}\t{zOTU_dict[zotu][1]}\t" + "\t".join(zOTU_dict[zotu][0]) + "\n")


# Call it for all relevant files
for file in os.listdir('.'):
    if file.endswith('.zotus.fasta'):
        sample = file.replace('.zotus.fasta', '')
        counts_file = f"{sample}_zotu_table.txt"
        fasta_file = f"{sample}.zotus.fasta"
        output_file = f"{sample}_zotu_table_with_seq.txt"
        add_seq_to_zotu(counts_file, fasta_file, output_file)

os.system("mkdir raw_zotu && mv *_zotu_table.txt raw_zotu && mv *zotus.fasta raw_zotu")
os.system("mkdir zotu_tables_with_sequences && mv *zotu_table_with_seq.txt zotu_tables_with_sequences")

###Creating a dictionary with sequences as keys and nested dictionary as values. Nested dictionary uses names of libraries as keys and number of reads as values:
seq_dict = {}
def get_headings(file):
   global headings
   headings = file.readline().strip('\t\n').split('\t')
   
def get_lib_name(heading):
    global lib_name
    for i in range(0, len(headings)):
       lib_name = headings[2]

def create_library(Lib_info):
   global lib_dict
   for line in Lib_info:
        LINE = line.strip('\t\n').split('\t')
        key = LINE[1]
        counts = LINE[-1]
        if not key in seq_dict.keys():
        #Creating nested dictionary with library names as keys and counts as values
            seq_dict[key] = {lib_name : counts}
        else:
            seq_dict[key].update({lib_name : counts})

# Store the current working directory
original_dir = os.getcwd()

# Change to the zotu_tables_with_sequences directory
zotu_tables_path = os.path.join(original_dir, "zotu_tables_with_sequences")
os.chdir(zotu_tables_path)

# Process files in the zotu_tables_with_sequences directory
files = os.listdir() 
for filename in files:
   with open(filename, "r") as Lib_info:
       get_headings(Lib_info)
       get_lib_name(headings)
       create_library(Lib_info)

# Change back to the original directory
os.chdir(original_dir)

# Creating a list with the names of the libraries
libs = []
for k1 in seq_dict.keys():
    libs += [lib for lib in seq_dict[k1] if lib not in libs] #List comprehension
    
### Creating an empty list which will be our final table
data = []
index = 0 #index of a table in data, basically 0 will raw with the first key, 1 - with the second ect
for seq in seq_dict.keys():
    data.append(["",seq, 0]) ### append  empty zotu_ID, sequence, total (which is zero at the beginning)
    for lib in libs: #For every library in our list
        data[index].append(0 if lib not in seq_dict[seq].keys() else seq_dict[seq][lib]) #append 0 if library is not in keys of nested dictionary for given sequence
        if lib in seq_dict[seq].keys(): ### summing Total for every library in the keys of sequence
            data[index][2] += int(seq_dict[seq][lib])
    index += 1 #increase index of one and go to the next sequence
 
def byTotal(Total): ###  Function which will indicate what to sort by
    return Total[2]
 
headers = ["#OTU_ID"] + libs ### Creating a list with headers
with open("all_libraries_zotu_table.txt", 'w') as bigFile:
    element = ""
    for header in headers:
        element += header + "\t" 
    bigFile.write(element[:-1] + '\n')
 
    data.sort(key=byTotal, reverse = True) ### sorting our table by Total with decreasing number of reads for each zOTU
    ###appending new zotu name for each sequence starting with the most abundant
    counter = 1
    for zotu in data:
        zotu[0] = "Zotu" + str(counter) 
        counter += 1
    for zotu in data:
        line = "" ###creating an empty string
        newData = [zotu[0]] + zotu[3:] ### adding ZOTU_ID and all the libraries to the list
        for element in newData:
            line += str(element) + "\t"
        bigFile.write(line[:-1] + '\n')

###Creating fasta file with sorted sequences of zOTUs:   
with open ("zotus.fasta", 'w') as fasta:
    for zotu in data:
        seq = ">" + zotu[0] + ";size=" + str(zotu[2]) + '\n' + zotu[1] + '\n'
        fasta.write(seq)
print("OK!")

# Run usearch clustering
subprocess.run("usearch -cluster_otus zotus.fasta -otus otus.fasta -relabel OTU -uparseout zotu_otu_relationships.txt -threads 20", shell=True)

# Run sed command to clean up zotus.fasta
subprocess.run("sed -E 's/;size=[0-9]+.*//g' zotus.fasta > new_zotus.fasta", shell=True)

os.system('clear')

# Run vsearch for taxonomic classification
subprocess.run(f"vsearch --sintax new_zotus.fasta -db {db_path} -tabbedout zotus.tax -strand both -sintax_cutoff 0.8 -threads {j}", shell=True)
subprocess.run(f"vsearch --sintax otus.fasta -db {db_path} -tabbedout otus.tax -strand both -sintax_cutoff 0.8 -threads {j}", shell=True)

# Remove redundant taxonomic labels
subprocess.run("sed -i 's/[dpcofgs]://g' zotus.tax", shell=True)
subprocess.run("sed -i 's/[dpcofgs]://g' otus.tax", shell=True)

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
with open("OTU_table.txt", "w") as bigFile:
    for LINE in data:
        print("\t".join(map(str, LINE)), file=bigFile)

# Run sed command to remove text inside parentheses
subprocess.run(["sed", "-i", "-E", "s/\\([^()]*\\)//g", "OTU_table.txt"])
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
# After processing, display the completion time
completion_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
print(f"Process completed at: {completion_time}")
