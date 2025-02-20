#usage: python3.7 Merging_sequencing_v2.py 000F-original_file_name1.seq
#001F-original_file_name2.seq 002R-original_file_name3.seq 003R-original_file_name4.seq

"""
Notice 1: First, all the Sanger sequence files need to be renamed and changed to plain
sequence format (seq format) with TextEdit by users.

Notice 2: The file names of the sequences to be merged starts with '000' plus 'F'
(forward sequencing result file) or 'R' (reverse sequencing result file),
and are in .seq format. This can be done by adding the four characters to the original
file name provided by the sequencing company.
The files are arranged in tandemly according to their real positions corresponding
to the sequencing plasmid or linear DNA.

Notice 3: The script requires Python 3.7.3, Biopython 1.7.4 and EMBOSS 6.6.0., please
install the three packages in advance before running the script.

Notice 4: This script is used to merge multiple successive Sanger DNA sequencing results.
It can merge 2 to hundreds of tandemly arranged Sanger sequencing results.
Please notice that a space key between each input item is obligatory, and the tab key can be
used to automatically fill up the input file name to simplify the input process. Also, the
backward slash "\" can be used for a multiline command input like this:

python3.7 Merge_Sanger_v2.py 001F_sequence_to_be_merged.seq \
002F_sequence_to_be_merged.seq 003F_sequence_to_be_merged.seq \
004R_sequence_to_be_merged.seq 005R_sequence_to_be_merged.seq \
006R_sequence_to_be_merged.seq
"""

import os,shutil
import sys
import string
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



for i in range(32):
    if i % 2 == 1:
        print(i*"(^_^)")
print("\n\n")


#get current path
current_path = os.getcwd()
print("Current path is %s." % current_path)
print("\n")

#make a new directory
folder_name = "merged_sequence"
dirs = os.listdir(current_path)
if folder_name not in dirs:
    os.mkdir(folder_name)
else:
    shutil.rmtree(folder_name)
    os.mkdir(folder_name)

number_of_files = len(sys.argv) -1
if len(sys.argv) -1 < 2:
    print("Must have 2 or more sequences. Please reinput sequences files to be merged.")
    sys.exit()
else:
    print("There are %s sequences to be joined." % number_of_files)

#copy sequence files to the new directory and change to F if it is in R
"""
F stands for the Sanger sequencing in forwarding direction and R for the Sanger sequencing
in reversing direction. It is convenient for the script to merge the Sanger sequences when
all the sequences have been changed to the forward sequencing direction. 
"""
for file in sys.argv[1:]:
    if file[3] is 'F':
        print(file)
        file_F = open(file)
        file_sequence = file_F.read()
        file_F.close()
        DNA_sequence_tmp_list = []
        
        for i in file_sequence:
            if i.isalpha():
                i = i.upper()
                DNA_sequence_tmp_list.append(i)
        DNA_sequence_tmp_str =''.join(DNA_sequence_tmp_list)        
        dir_sub = os.path.join(current_path,folder_name)
        os.chdir(dir_sub)
        f = open(file[0:4] + '.seq','w')
        f.write(DNA_sequence_tmp_str)
        f.close()
        os.chdir(current_path)
        
    if file[3] is 'R':
        print(file)
        file_R = open(file)
        file_sequence = file_R.read()
        file_R.close()
        DNA_sequence_tmp_list = []
        
        for i in file_sequence:
            if i.isalpha():
                i = i.upper()
                DNA_sequence_tmp_list.append(i)
        DNA_sequence_tmp_str =''.join(DNA_sequence_tmp_list)
        DNA_sequence_tmp_Seq = Seq(DNA_sequence_tmp_str)
        DNA_sequence_tmp_Seq_F = DNA_sequence_tmp_Seq.reverse_complement()
        DNA_sequence_tmp_Seq_F_str = str(DNA_sequence_tmp_Seq_F)
        
        dir_sub = os.path.join(current_path,folder_name)
        os.chdir(dir_sub)
        f = open(file[0:4] + '.seq','w')
        f.write(DNA_sequence_tmp_Seq_F_str)
        f.close()
        os.chdir(current_path)
    

#function for the boundaries
"""
The function below parses the output of the EMBOSS needle alignment between the tandemly
arranged Sanger sequencing files to obtain the positions of the first 50 nucleotides match.
The EMBOSS needle shows the alignment result in segments of 50 nucleotides.
This boundary is used to determine which part of each Sanger sequence is kept for joining
to the final merged sequence. This operation is actually the same as that when we merge the
Sanger sequencing files manually.
"""
def needle_align_to_get_boundaries(f1,f2):
    alignment_a_left = 0
    alignment_b_left = 0
    output_file_name = f1.split('.')[0] + f2.split('.')[0] + ".needle"
    
    needle_cline = NeedleCommandline()
    needle_cline.asequence = f1
    needle_cline.bsequence = f2
    needle_cline.gapopen = 10
    needle_cline.gapextend= 0.5
    needle_cline.outfile = output_file_name
    print(needle_cline)
    stdout, stderr = needle_cline()
    print(stdout + stderr)

    #open the needle alignment output file and get boundaries
    file = open(output_file_name)
    file_lines = file.readlines()
    file.close()

    for line in file_lines:
        print(line, end="")

    alignment_a_squence_positions = []
    alignment_b_squence_positions = []

    file = open(output_file_name)

    new_line1 = file.readline()
    new_line2 = file.readline()

    while len(new_line2):
        line_a = new_line1
        line_b = new_line2
        new_line2 = new_line2.strip()
        
        if (50*'|' in new_line2):

            line_b = file.readline()
      
            alignment_a_squence_line_str = line_a.strip()
            alignment_b_squence_line_str = line_b.strip()
            print("The beginning of excellent alignment is shown below.\n")
            
            alignment_a_squence_line_str_split = alignment_a_squence_line_str.split()
            print(alignment_a_squence_line_str_split[0].ljust(5,' '),\
                  alignment_a_squence_line_str_split[1],\
                  alignment_a_squence_line_str_split[2].rjust(6,' '),\
                  sep="")
            alignment_b_squence_line_str_split = alignment_b_squence_line_str.split()
            print(alignment_b_squence_line_str_split[0].ljust(5,' '),\
                  alignment_b_squence_line_str_split[1],\
                  alignment_b_squence_line_str_split[2].rjust(6,' '),\
                  sep="")
     
            print("\n")

            alignment_a_left = int(alignment_a_squence_line_str.split()[0])
            alignment_b_left = int(alignment_b_squence_line_str.split()[0])
            break
        else:
            new_line1 = new_line2         #notice the skill here, we must go step by step through the lines
            new_line2 = file.readline()
    file.close()

    return alignment_a_left, alignment_b_left

#end of function

#align with needle progressively to get boundaries
dir_sub = os.path.join(current_path,folder_name)
os.chdir(dir_sub)
dir_new = os.listdir()
dir_new.sort()

Es_list = []
Es_list.append(1)

for i in range(len(dir_new) - 1):
    E = needle_align_to_get_boundaries(dir_new[i],dir_new[i+1])
    Es_list.append(E[0])
    Es_list.append(E[1])

#length of the last sequence
f1 = open(dir_new[-1])
file_sequence = f1.read()
f1.close()
DNA_sequence_temp_list = []
for i in file_sequence:
    if i.isalpha():
        i = i.upper()
        DNA_sequence_tmp_list.append(i)
length_of_last_file_sequence = len(DNA_sequence_tmp_list)

Es_list.append(length_of_last_file_sequence + 1)

#till now, got the boundaries for each sequence

#list to store the final merged sequence
merged_sequence_list = []        
        
for i in range(len(dir_new)):
    f = open(dir_new[i])
    file_sequence = f.read()
    f.close()
    DNA_sequence_tmp_list = []
    for j in file_sequence:
                if j.isalpha():
                    j = j.upper()
                    DNA_sequence_tmp_list.append(j)
    merged_sequence_list += DNA_sequence_tmp_list[Es_list[2*i]-1:Es_list[2*i+1]-1]
"""
The content of Es_list, 2n elements in total, where n=len(dir_new) is the number of
sequence files to be merged.

E0=1,E1;    E2,E3;    E4,E5;    E6,E7;    ...,    E2n-2,E2n-1
S0          S1        S2        S3        ...,    Sn-1

The above S0,S1,...,Sn-1 denote the sequence files to be merged.
The E0,E1,E2,...,E2n-1, which are stored in list Es_list, represent the coordinates
of intervals to join to full-length sequence. The alignment program needle
starts nucleotide numbers from 1.
However, in Python the list starts from 0. Thus, a mathematical interval [X,Y), which is
equal to a mathematical interval [X,Y-1], corresponds to the list elements from
X-1 to Y-2, i.e., list[X-1:Y-1].
Thus, DNA_sequence_tmp_list[Es_list[2*i]-1:Es_list[2*i+1]-1] stands for the nucleotides
that will be used for merging for each Sanger sequence file that has been changed to
forward sequencing.
"""

#print the merged sequence
merged_sequence_str = ''.join(merged_sequence_list)
print("The merged DNA sequence is shown below.\n")
for j in range(0,len(merged_sequence_str),100):
    DNA_string_100_per_line_str = merged_sequence_str[j:j+100]
    print(DNA_string_100_per_line_str)
print("\n")

#write to the file
merged_file_name = "merged" + ".seq"
file=open(merged_file_name,'w')
file.write(merged_sequence_str)
file.close()
os.system('rm *.needle')
os.system('rm 00*')


