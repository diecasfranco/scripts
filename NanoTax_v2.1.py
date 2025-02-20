#! /usr/bin/env python3
# Authors: Diego Franco, Junchen Deng
# Affiliation: Symbio Group, Institute of Enviromental Sciences, Jagiellonian University, Poland
# Version 2.1

# requirement: blast/diamond, samtools, minimap2, bowtie2, pileup.sh

import os, argparse, itertools, sys, csv, time

parser = argparse.ArgumentParser(
			description='This script is intended to produce a table with both contig information (e.g. average coverage, GC content) and the corresponding taxonomy for output contigs from long-reads and short-reads assemblers (e.g. Canu, Flye for Nanopore, Megahit for illumina)'	
)

# positional arguments
parser.add_argument("contigs", metavar='<contigs>', help="the path to the contig fasta file/folder")
parser.add_argument("database", metavar='<database.FASTA>', type=str, help="the path to the database")
parser.add_argument("db_type", metavar='<database type>', type=str, help="database type (option: 'protein' 'dna')")

# optional arguments
parser.add_argument("-db_2", "--database_2", metavar='<.FASTA>', type=str, help="the path to the 2nd database. The db_type of the 2nd database has to be different from the 1st database")
parser.add_argument("-bt", "--blast_tool", metavar='<blast/diamond>', type=str, help="the type of tools for blast (default: %(default)s) (option: 'blast', 'diamond')", default="blast")
parser.add_argument("-blastx_task", "--blastx_task", metavar='<str>', type=str, help="the task type of blastx (option: 'blastx' 'blastx-fast') (default: %(default)s)", default="blastx")
parser.add_argument("-blastn_task", "--blastn_task", metavar='<str>', type=str, help="the task type of blastn (option: 'blastn' 'megablast') (default: %(default)s)", default="megablast")
parser.add_argument("-diamond_sen", "--diamond_sensitivity", metavar='<str>', type=str, help="the sensitivity of diamond (option: 'fast' 'mid-sensitive' 'sensitive' 'more-sensitive' 'very-sensitive' 'ultra-sensitive') (default: %(default)s)", default="sensitive")
parser.add_argument('-BAM', '--BAM', metavar='', help="path to BAM file; enable skipping the mapping step")
parser.add_argument('-SAM', '--SAM', metavar='', help="path to SAM file; enable skipping the mapping step")
parser.add_argument('-cov', '--coverage', metavar='<from_pileup.sh>', help="path to coverage fiel from pileup.sh; enable to skip both mapping and pileup.sh")
parser.add_argument('-bnf', '--blastn_file', metavar='<outfmt6>', help="path to blastn output; enable skipping the blastn step")
parser.add_argument('-bxf', '--blastx_file', metavar='<outfmt6>', help="path to blastx output; enable skipping the blastx step")
parser.add_argument('--truecontigs', action='store_true', help="enable analysing only contigs with blast hits")
parser.add_argument("-ONT_fastq", "--ONT_fastq", metavar='<FASTQ>', help="the path to the Nanopore reads fastq file/folder")
parser.add_argument("-r1", "--r1", metavar='<FASTQ>', help="the path to the r1 pair-end reads file/folder")
parser.add_argument("-r2", "--r2", metavar='<FASTQ>', help="the path to the r2 pair-end reads file/folder")
parser.add_argument("-o", "--output_dir", metavar='', type=str, help="output directory name (default: %(default)s)", default="Assigned_Taxonomy")
parser.add_argument("-prefix", "--prefix", metavar='', type=str, help="prefix of each output file name (default: %(default)s)", default="output")
parser.add_argument('-c', '--cores', metavar='', type=int,
			help='The number of CPU cores the script will use (default: %(default)s)', default=8)

# if no arguments were given, printing the help message (args = "--help")
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

contigs = args.contigs
db = args.database
db_type = args.db_type

db_2 = args.database_2
bt = args.blast_tool
blastx_task = args.blastx_task
blastn_task = args.blastn_task
diamond_sen = args.diamond_sensitivity
ONT_fastq = args.ONT_fastq
BAM = args.BAM
SAM = args.SAM
coverage = args.coverage
blastn_file = args.blastn_file
blastx_file = args.blastx_file
truecontigs = args.truecontigs
r1_fastq = args.r1
r2_fastq = args.r2
outdir = args.output_dir
prefix = args.prefix
cores = args.cores

output = "{0}/{1}".format(outdir,prefix) 

start_time = time.time()
os.system("mkdir {0}".format(outdir))    # default = ./Assigned_Taxonomy"

def ImportTax(path):
    table_tax = {}
    with open(path, "r") as Table:
        for line in Table: 
            LINE = line.strip().split()
            header = LINE[0]
            taxonomy = LINE[1]
            table_tax[header] = taxonomy
    return(table_tax)

# 1) assign taxonomy to each contig by blast or diamond
# output: prefix.blastn and/or prefix.blastx
if db_type == 'dna':
    print("\n\nThe input db is DNA")
    print(">Performing blastn...")
    if blastn_file != None:
        print("blastn file provided; skipping blastn")
        os.system("cp {0} {1}.blastn".format(blastn_file, output))
    else:
        os.system("blastn -task {4} -query {0} -db {1} -outfmt 6 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads {2} > {3}.blastn".format(contigs, db, cores, output, blastn_task))
        print("blastn done\n")
else:
    print("\n\nThe input db is protein")
    if blastx_file != None:
        print("blastx file provided; skipping blastx")
        os.system("cp {0} {1}.blastx".format(blastx_file, output))
    else:
        if bt == "diamond":
            print(">Performing diamond blastx... ")
            os.system("diamond blastx -p {2} --db {1}.dmnd --query {0} -o {3}.blastx --evalue 1e-10 --outfmt 6 --max-target-seqs 1 --{4}".format(contigs, db.split(".fasta")[0], cores, output, diamond_sen))
            print("diamond blastx done\n")
        else:
            print(">Performing blastx... ")
            os.system("blastx -task {4} -query {0} -db {1} -outfmt 6 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads {2} > {3}.blastx".format(contigs, db, cores, output, blastx_task))
            print("blastx done\n")
     
if db_2 != None:
    if db_type == 'protein':
        print("\n\nThe 2nd input db is DNA")
        print(">Performing blastn...")
        if blastn_file != None:
            print("blastn file provided; skipping blastn")
            os.system("cp {0} {1}.blastn".format(blastn_file, output))
        else:
            os.system("blastn -task {4} -query {0} -db {1} -outfmt 6 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads {2} > {3}.blastn".format(contigs, db_2, cores, output, blastn_task))
            print("blastn done\n")
    else:
        print("\n\nThe 2nd input db is protein")
        if blastx_file != None:
            print("blastx file provided; skipping blastx")
            os.system("cp {0} {1}.blastx".format(blastx_file, output))
        else:
            if bt == "diamond":
                print(">Performing diamond blastx... ")
                os.system("diamond blastx -p {2} --db {1}.dmnd --query {0} -o {3}.blastx --evalue 1e-10 --outfmt 6 --max-target-seqs 1 --{4}".format(contigs, db_2.split(".fasta")[0], cores, output, diamond_sen))
                print("diamond blastx done\n")
            else:
                print(">Performing blastx... ")
                os.system("blastx -task {4} -query {0} -db {1} -outfmt 6 -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads {2} > {3}.blastx".format(contigs, db_2, cores, output, blastx_task))
                print("blastx done\n")
# make sure the provided blast files are in the output directory
if blastn_file != None:
    os.system("cp {0} {1}.blastn".format(blastn_file, output))
if blastx_file != None:
    os.system("cp {0} {1}.blastx".format(blastx_file, output))

# 2) import taxonomy table as a dictionary
blastn_count = 0
blastx_count = 0
if os.path.exists("{}.blastn".format(output)): 
    blastn_count = 1
    table_blastn = ImportTax("{}.blastn".format(output))

if os.path.exists("{}.blastx".format(output)):
    blastx_count = 1
    table_blastx = ImportTax("{}.blastx".format(output))

#####################################
if coverage != None:
    print("coverage file provided; skip mapping")
    os.system("cp {0} {1}.coverage".format(coverage, output))
else:
    if SAM != None:
        print("SAM file provided; skip mapping, produce coverage file, and convert SAM to BAM")
        os.system("pileup.sh in={0} out={1}.coverage".format(SAM, output))
        print("Coverage done\n")
        os.system("samtools view -b -@ 4 {0} | samtools sort -@ 4 - -o {1}.bam".format(SAM, output))
        os.system("samtools index {0}.bam".format(output))
    else:
        if BAM != None:
            print("BAM file provided; skip mapping and produce coverage file")
            os.system ("pileup.sh in={0} out={1}.coverage".format(BAM, output))
            print("Coverage done\n")
        else:
            ## 3) [enabled when "--truecontigs" is set] output contigs with blast hit
            if truecontigs: 
                with open(contigs, "r") as FASTA:
                    contig_dict = {}
                    contig = ""
                    for line in FASTA:
                        if line.startswith(">"):
                            if contig != "":    # save the existing header and sequence before overwritting 
                                contig_dict[header] = contig
                                contig = ""
                            header = line.strip("\n").split()[0][1:]    # remove the ">" from the header for the following steps
                        else:
                            contig += line.strip("\n").upper()
                    contig_dict[header] = contig

                with open("{}.TrueContigs.fasta".format(output), "w") as OUTPUT:
                    for header in contig_dict.keys():
                        if blastn_count == 1 and blastx_count == 1:
                            if header in table_blastn.keys() or header in table_blastx.keys():
                                print(">"+ header, "\n", contig_dict[header], "\n", end = "", sep = "", file = OUTPUT)
                        elif blastn_count == 1:
                            if header in table_blastn.keys():
                                print(">"+ header, "\n", contig_dict[header], "\n", end = "", sep = "", file = OUTPUT)
                        elif header in table_blastx.keys():
                            print(">"+ header, "\n", contig_dict[header], "\n", end = "", sep = "", file = OUTPUT)

            ## 4) mapping reads
            if ONT_fastq != None: 
                print("Mapping long reads with Minimap2...")
                if truecontigs: 
                    os.system("minimap2 -ax map-ont {2}.TrueContigs.fasta {0} -t {1} > {2}.sam".format(ONT_fastq, cores, output))
                else:
                    os.system("minimap2 -ax map-ont {3} {0} -t {1} > {2}.sam".format(ONT_fastq, cores, output, contigs))   
                print("Minimap2 done\n")
                
            if r1_fastq != None and r2_fastq != None:
                print("Mapping pair-end short reads with bowtie2...")
                if truecontigs:
                    os.system("bowtie2-build {0}.TrueContigs.fasta {0}".format(output))
                    os.system("bowtie2 --threads {2} -x {3} -1 {0} -2 {1} --no-unal | samtools sort -@ 4 - -o {3}.bam".format(r1_fastq, r2_fastq, cores, output))
                else:
                    os.system("bowtie2-build {0} {1}".format(contigs, output))
                    os.system("bowtie2 --threads {2} -x {3} -1 {0} -2 {1} --no-unal | samtools sort -@ 4 - -o {3}.bam".format(r1_fastq, r2_fastq, cores, output))
                print("bowtie2 done\n")    
                
            os.system ("pileup.sh in={0}.bam out={0}.coverage".format(output))
            print("Coverage done\n")
            os.system("samtools index {0}.bam".format(output))
            os.system("rm {0}/*.bt2".format(outdir))

## 5) combining information in taxonomy and coverage tables
# the first line will be categories, i.e. "#ID, Avg_fold, Length, Ref_GC, Covered_percent, Covered_bases, Plus_reads, Minus_reads, Read_GC, Median_fold, Std_Dev"
with open("{}.coverage".format(output), "r") as Table:
    Table_cov = []
    for line in Table: 
        LINE = line.strip().split()
        Table_cov.append(LINE)
row_count = len(Table_cov)

if blastn_count == 1: 
    Table_cov[0].append("blastn")
    for row in range(1, row_count):
        if Table_cov[row][0] in table_blastn.keys(): 
            header = Table_cov[row][0]
            taxonomy = table_blastn[header]
            Table_cov[row].append(taxonomy)
        else: 
            Table_cov[row].append("NA")

if blastx_count == 1: 
    Table_cov[0].append("blastx")
    for row in range(1, row_count):
        if Table_cov[row][0] in table_blastx.keys(): 
            header = Table_cov[row][0]
            taxonomy = table_blastx[header]
            Table_cov[row].append(taxonomy)
        else: 
            Table_cov[row].append("NA")

# provide a summerized taxonomy       
Table_cov[0].append("Taxonomy")
if blastn_count + blastx_count == 1:
    for row in range(1, row_count):
        tax = Table_cov[row][-1].split("_")[0]    # only take the header before "_" as the taxonomy
        Table_cov[row].append(tax)
            
        if "mito" in Table_cov[row][-2]:
            Table_cov[row][-1] = "mitochondrion"
        if Table_cov[row][-1] == "Acyrthosiphon": 
            Table_cov[row][-1] = "Host"
            
if blastn_count == 1 and blastx_count == 1:    
    for row in range(1, row_count):
        tax1 = Table_cov[row][-1].split("_")[0]
        tax2 = Table_cov[row][-2].split("_")[0]
        if tax1 == tax2: 
            Table_cov[row].append(tax1)
        if tax1 == 'NA' and tax2 != 'NA':
            Table_cov[row].append(tax2)
        if tax2 == 'NA' and tax1 != 'NA':
            Table_cov[row].append(tax1)
        if tax1 != 'NA' and tax2 != 'NA' and tax1 != tax2:
            Table_cov[row].append("conflicting")
        
        if "mito" in Table_cov[row][-3] or "mito" in Table_cov[row][-2]:
            Table_cov[row][-1] = "mitochondrion"
        if Table_cov[row][-1] == "Acyrthosiphon": 
            Table_cov[row][-1] = "Host"
                      
# output table
if blastn_count + blastx_count == 1:
    if blastn_count == 1: 
        col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "blastn"]
    else:
        col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "blastx"]

if blastn_count == 1 and blastx_count == 1: 
    col_keep = ["#ID", "Avg_fold", "Length", "Read_GC", "blastn", "blastx", "Taxonomy"]
    
col_pos = []
for name in col_keep: 
    col_pos.append(Table_cov[0].index(name))

with open("{}.taxonomy".format(output), "w") as OUTPUT:
    for row in Table_cov: 
        col_position = 0
        for col in row[:-1]:
            if col_position in col_pos:
                print(col, file = OUTPUT, end = "\t")
            col_position += 1
        print(row[-1], file = OUTPUT)

end_time = time.time()
print("running time: {} seconds".format(end_time - start_time))     
print("Na zdrowie! Salud! 乾杯! gānbēi(干杯^-^)! Saude! Cheers! Skål!")



