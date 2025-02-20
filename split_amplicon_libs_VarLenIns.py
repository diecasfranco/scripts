#! /usr/bin/env python3

import sys, os, re
from Bio.Data import IUPACData
from Bio import SeqIO

if len(sys.argv) != 4:
	sys.exit('\nsplit_amplicon_libs_VarLenIns.py v. 0.2. Piotr Lukasik, 4th August 2020\n'
	         '-------------------------------------------------------------\n'
	         'This script splits up a pair of fastq files corresponding to a multiplex amplicon library '
	         'into a collection of fastq file pairs corresponding to different targeted regions,\n'
	         'using the sequences of primers for alternative targets provided. The script considers'
	         'that the primer sequences may be preceded by up to four other bases (variable-length insert).\n'
	         '    Usage: ./split_amplicon_libs_VarLenIns.py <list_of_targets> <list_of_libraries> <output_dir> \n'
	         '    e.g., ./split_amplicon_libs_VarLenIns.py targets.txt samples.txt ~/test_data')
	         
Script, Target_List, Sample_List, Output_Dir = sys.argv

print('\nsplit_amplicon_libs_VarLenIns.py v. 0.2. Piotr Lukasik, 4th August 2020\n'
	   '-------------------------------------------------------------\n')


###################
### Block 1. Reading and converting amplicon target list ....
###################

print("Reading the list of targets .....", end="")

Target_List_Imported = []
TARGET_LIST = open(Target_List, 'r', encoding="utf-8")
for line in TARGET_LIST:
    if not line.startswith("#") and len(line) > 0:
        if len(line.strip().split()) != 3:
            print("There is a problem with this line, which will be skipped:\n    ", end = "")
            print(line.strip())
        else:
            Target_List_Imported.append(line.strip().split())
TARGET_LIST.close()
# Target_List_Imported = [[target_name, primer_F, primer_R], [target_name2, primer_F, primer_R]...]
# contains lists of targeted amplicon regions: target name, forward primer, reverse primer

		
Target_List = []
### Now, I am going to convert primers (with ambiguities) to search terms for re.match(), such as 'AG[AG]CT[CTG]A'
for line in Target_List_Imported:
    if len(line) == 0:
        break
    #primer_R1_ambig = ""      # OK if we provide primers without any variable length inserts, etc.
    primer_R1_ambig = ".{0,4}" # assumes that there can be up to four of whatever prior to the sequence!
    #for nt in line[1]: 
    for nt in line[1][1:]:     # skips the first base in the primer!
        value = IUPACData.ambiguous_dna_values[nt] 
        if len(value) == 1: 
            primer_R1_ambig += value 
        else: 
            primer_R1_ambig += "[%s]" % value 
    #primer_R2_ambig = ""        # OK if we provide primers without any variable length inserts, etc.
    primer_R2_ambig = ".{0,4}"   # assumes that there can be up to four of whatever prior to the sequence!
    #for nt in line[2]: 
    for nt in line[2][1:]:       # skips the first base in the primer!
        value = IUPACData.ambiguous_dna_values[nt] 
        if len(value) == 1: 
            primer_R2_ambig += value 
        else: 
            primer_R2_ambig += "[%s]" % value 
                
    Target_List.append([line[0], primer_R1_ambig, primer_R2_ambig, []])
    # target_name, primer_F, primer_R, list of sequences classified as corresponding to the target. Empty so far.

Target_List.append(["Unclassified", '', '', []])

print("OK! %s targets provided:\n    " % str(len(Target_List)-1), end = "")
for target_info in Target_List[:-1]:
    print(target_info[0], end = ", ")
print("")



###################
### Block 2. Reading sample list ....
###################

print("\nReading the list of samples .....", end ="")

Sample_List_Imported = []
SAMPLE_LIST = open(Sample_List, 'r', encoding="utf-8")
for line in SAMPLE_LIST:
    if not line.startswith("#") and len(line) > 0:
        if len(line.strip().split()) != 3:
            print("There is a problem with this line, which will be skipped:\n    ", end = "")
            print(line.strip())
        else:
            Sample_List_Imported.append(line.strip().split())
SAMPLE_LIST.close()
# Sample_List_Imported = [[libr1, libr1_R1.fastq, libr1_R2.fastq], ...]
# contains lists of targeted amplicon regions: target name, forward primer, reverse primer


Target_Counts_per_Sample = [['']]
for target_info in Target_List:
    Target_Counts_per_Sample[0].append(target_info[0])
# This list is going to be used for keeping summary data for libraries

print("OK! %d samples will be sorted" % len(Sample_List_Imported))



###################
### Block 3. Processing samples ....
###################


for library in Sample_List_Imported:
    print("Sorting data for sample %s ................." % library[0], end = "")
    
    for target_info in Target_List:
        target_info[3] = []
    
    
    ### Now, reading in R1 and R2 fastq files. 
    READ1 = SeqIO.to_dict(SeqIO.parse(library[1], "fastq"))
    READ2 = SeqIO.to_dict(SeqIO.parse(library[2], "fastq"))
    
    ### Now, creating a 'Read_db' containing {Seq_ID : [R1_sequence, R2_sequence]}
    Read_db = {}
    for key in READ1:
        Read_db[key] = [str(READ1[key].seq)]
    for key in READ2:
        Read_db[key].append(str(READ2[key].seq))
    

    ### Now, comparing target primers with beginnings of sequences, and adding the identified Seq_IDs to lists associated with correct targets
    for key, value in Read_db.items():   # for every read, e.g., {'Seq1234' : ['ACGTGTTG', 'CTTGTGAT']}:
        for target_info in Target_List:  # for every target, e.g., ['16Sv4', 'CTYC', 'CAAC']
            if re.match(target_info[1], value[0]) and re.match(target_info[2], value[1]): # if the primer seqs match the beginnings of target seqs,
                target_info[3].append(key)   # add sequence ID to the list of sequences matching the target
                break

    ### Now, exporting reads, while adding info to 'Target_Counts_per_Sample'... 
    for target_info in Target_List:
        R1_export_fastq = open(Output_Dir + "/" + target_info[0] + "_" + library[0] + "_R1.fastq", "w")
        R2_export_fastq = open(Output_Dir + "/" + target_info[0] + "_" + library[0] + "_R2.fastq", "w")
        for ID in target_info[3]:
            SeqIO.write(READ1[ID], R1_export_fastq, "fastq")
            SeqIO.write(READ2[ID], R2_export_fastq, "fastq")
        
        R1_export_fastq.close()
        R2_export_fastq.close()
    
    Current_Target_Counts = [library[0]]
    for target_info in Target_List:
        Current_Target_Counts.append(len(target_info[3]))
    Target_Counts_per_Sample.append(Current_Target_Counts)


    print("OK! %d paired-end sequences classified!" % len(Read_db))



###################
### Block 4. Almost done, just printing summaries :)
###################

Summary_file = open(Output_Dir + "/000_splitting_summary.txt", "w")

print("\nJob complete!\n")
print("Target counts in different samples:")
for line in Target_Counts_per_Sample:
    for item in line:
        print(item, end = "\t")
        print(item, end = "\t", file = Summary_file)
    print("")
    print("", file = Summary_file)

            



