# -*- coding: utf-8 -*-
"""
Python 3.10

@author: Andy Vierstraete

Trim adapters from Oxford Nanopore sequencing reads.
Demultiplex reads based on default or custom barcodes.
De-UMI-plex reads based on used UMI barcodes.

depends on:
Lightweight, super fast C/C++ library for sequence alignment using edit (Levenshtein) 
distance:
    https://pypi.org/project/edlib/#description
    python3 -m pip install edlib

"""
import edlib
import time
from Bio import SeqIO
import sys
import os
from threading import Thread, Event
from multiprocessing import Process, Lock, Queue
from multiprocessing.sharedctypes import Value
import argparse
from itertools import zip_longest
import glob
import csv
import pickle
import re
import urllib.request
import gzip
#==============================================================================
"""
Abbreviations from nanopore:
    
LIGATION KITS

  lsk: ligation sequencing kit (lsk109, lsk110, lsk114)
  nbd: native barcoding kit           (NB barcodes)
  mlk: multiplex ligation kit         (NB barcodes)
  cs9: Cas9 sequencing kit                              Discontinued
  
  pbc: PCR barcoding expansion kit    (BC barcodes) (not a sequencing kit)
  
RAPID KITS
  
  rad: rapid sequencing kit
  rbk: rapid barcoding kit            (RB barcodes)
  ulk: ultra long DNA sequencing kit
  lrk: Field sequencing kit                             Discontinued
  psk: PCR sequencing kit                               Discontinued
  pbk: PCR barcoding kit              (BP barcodes)     Discontinued
  rpb: rapid PCR barcoding kit        (RLB barcodes)
  16s: 16s barcoding kit              (16s barcodes)
  rab: 16s rapid barcoding kit        (16s barcodes)    Discontinued

RNA - cDNA KITS

  rna: direct RNA sequencing kit
  pcs: cDNA PCR sequencing kit   (contains UMI ???)
  dcs: direct cDNA sequencing kit                       Discontinued
  pcb: cDNA PCR barcoding kit         (BP barcodes)

"""
adapters = {
    'ATCGCCTACCGTGAC---BC---AGAGTTTGATCMTGGCTCAG': ['flank_16s_rab_forward'], 
    'ATCGCCTACCGTGAC---BC---CGGTTACCTTGTTACGACTT': ['flank_16s_rab_reverse'], 
    'GGTGCTG---BC---TTAACCT': ['flank_pbc_forward'], 
    'GCTTGGGTGTTTAACC---BC---GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA': ['flank_rbk'], 
    'ATCGCCTACCGTGAC---BC---ACTTGCCTGTCGCTCTATCTTC': ['flank_pbk_forward'], 
    'ATCGCCTACCGTGAC---BC---TTTCTGTTGGTGCTGATATTGC': ['flank_pbk_reverse'], 
    'ATCGCCTACCGTGA---BC---TTGCCTGTCGCTCTATCTTC': ['flank_pcb_forward'], 
    'ATCGCCTACCGTGA---BC---TCTGTTGGTGCTGATATTGC': ['flank_pcb_reverse'], 
    'ATCGCCTACCGTGAC---BC---CGTTTTTCGTGCGCCGCTTC': ['flank_rpb'], 
    'AAGGTTAA---BC---CAGCACCT': ['flank_nbd'], 
    'MMTGTACTTCGTTCAGTTACGTATTGCT': ['adapter_lsk109_front', 'adapter_lsk114_front', 'adapter_nbd_front',
                                     'adapter_rad_front', 'adapter_cs9_front', 'adapter_cdna_pcs_rapt_front', 
                                     'adapter_cdna_dcs_amx_front', 'adapter_mlk_front', 'adapter_rbk_front',
                                     'adapter_ulk_front', 'adapter_lrk_front', 'adapter_psk_front', 'adapter_pbk_front',
                                     'adapter_rpb_front', 'adapter_16s_front', 'adapter_rab_front', 'adapter_pcb_front'], 
    'AGCAATACGTAACTGAACGAAGT': ['adapter_lsk109_rear'], 
    'MMTGTACTTCGTTCAGTTACGTATTGC': ['adapter_lsk110_front'], 
    'AGCAATACGTAACTGAAC': ['adapter_lsk110_rear'], 
    'GCAATACGTAACTGAACGAAGTACAGG': ['adapter_lsk114_rear', 'adapter_mlk_rear'], 
    'ACGTAACTGAACGAAGTACAGG': ['adapter_nbd_rear'], 
    'ACTTCGTTCAGTTACGTATTGCT': ['adapter_middle'], 
    'GCAATACGTAACTGAACGAAGT': ['adapter_cs9_rear', 'adapter_cdna_dcs_amx_rear'], 
    'GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC': ['adapter_rna_rta_front'],  # reverse transcription adapter
    'GAGGCGAGCGGTCAATTTTCCTAAGAGCAAGAAGAAGCC': ['adapter_rna_rta_rear'], # reverse transcription adapter
    'ATGATGCAAGATACGCAC': ['adapter_rna_rmx_front'], # rna adapter mix
    'GAGGCGAGCGGTCAATTTGCAATATCAGCACCAACAGAAACAACCATCGTCTATCCCTCATCATCAGAACCTACTA': ['adapter_rna_rmx_rear'],  # rna adapter mix
    'CTTGCCTGTCGCTCTATCTTCAGAGGAG': ['adapter_cdna_pcs_rtp_front'], # reverse transcription primer
    'TTTCTGTTGGTGCTGATATTGCTTTVVVVTTVVVVTTVVVVTTVVVVTTTGGG': ['adapter_cdna_pcs_sspii_front', 
                                                              'adapter_cdna_dcs_ssp_front'], # strand switching primer
    'TTTCTGTTGGTGCTGATATTGCTGGG': ['adapter_cdna_dcs_ssp_front'], # strand switching primer
    'ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGACTTGCCTGTCGCTCTATCTTC': ['adapter_cdna_pcs_front'], # cDNA primer
    'ATCGCCTACCGTGACAAGAAAGTTGTCGGTGTCTTTGTGTTTCTGTTGGTGCTGATATTGC': ['adapter_cdna_pcs_rear'], # cDNA primer
    'CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG': ['adapter_cdna_pcs_rt_front'], # reverse transcription adapter
    'ACTTGCCTGTCGCTCTATCTTC': ['adapter_cdna_dcs_vnp_front'], 
    'CACAAAGACACCGACAACTTTCTT': ['bc01_nbd_mlk'],  # from ligation kits nbd
    'ACAGACGACTACAAACGGAATCGA': ['bc02_nbd_mlk'], 
    'CCTGGTAACTGGGACACAAGACTC': ['bc03_nbd_mlk'], 
    'TAGGGAAACACGATAGAATCCGAA': ['bc04_nbd_mlk'], 
    'AAGGTTACACAAACCCTGGACAAG': ['bc05_nbd_mlk'], 
    'GACTACTTTCTGCCTTTGCGAGAA': ['bc06_nbd_mlk'], 
    'AAGGATTCATTCCCACGGTAACAC': ['bc07_nbd_mlk'], 
    'ACGTAACTTGGTTTGTTCCCTGAA': ['bc08_nbd_mlk'], 
    'AACCAAGACTCGCTGTGCCTAGTT': ['bc09_nbd_mlk'], 
    'GAGAGGACAAAGGTTTCAACGCTT': ['bc10_nbd_mlk'], 
    'TCCATTCCCTCCGATAGATGAAAC': ['bc11_nbd_mlk'], # from rapid(RBK) PCR(PBC) PCR(PBK) rapid_PCR(RPB) 16S(RAB_16S) 
    'TCCGATTCTGCTTCTTTCTACCTG': ['bc12_nbd_mlk'], # rapid_V14(RBK) cDNA(PCB) kits
    'AGAACGACTTCCATACTCGTGTGA': ['bc13_nbd_mlk', 'bc13_rbk_pbc_rab_16s_pcb'],
    'AACGAGTCTCTTGGGACCCATAGA': ['bc14_nbd_mlk', 'bc14_rbk_pbc_rab_16s_pcb'], 
    'AGGTCTACCTCGCTAACACCACTG': ['bc15_nbd_mlk', 'bc15_rbk_pbc_rab_16s_pcb'], 
    'CGTCAACTGACAGTGGTTCGTACT': ['bc16_nbd_mlk', 'bc16_rbk_pbc_rab_16s_pcb'], 
    'ACCCTCCAGGAAAGTACCTCTGAT': ['bc17_nbd_mlk', 'bc17_rbk_pbc_rab_16s_pcb'], 
    'CCAAACCCAACAACCTAGATAGGC': ['bc18_nbd_mlk', 'bc18_rbk_pbc_rab_16s_pcb'], 
    'GTTCCTCGTGCAGTGTCAAGAGAT': ['bc19_nbd_mlk', 'bc19_rbk_pbc_rab_16s_pcb'], 
    'TTGCGTCCTGTTACGAGAACTCAT': ['bc20_nbd_mlk', 'bc20_rbk_pbc_rab_16s_pcb'], 
    'GAGCCTCTCATTGTCCGTTCTCTA': ['bc21_nbd_mlk', 'bc21_rbk_pbc_rab_16s_pcb'], 
    'ACCACTGCCATGTATCAAAGTACG': ['bc22_nbd_mlk', 'bc22_rbk_pbc_rab_16s_pcb'], 
    'CTTACTACCCAGTGAACCTCCTCG': ['bc23_nbd_mlk', 'bc23_rbk_pbc_rab_16s_pcb'], 
    'GCATAGTTCTGCATGATGGGTTAG': ['bc24_nbd_mlk', 'bc24_rbk_pbc_rab_16s_pcb'], 
    'GTAAGTTGGGTATGCAACGCAATG': ['bc25_nbd_mlk', 'bc25_rbk_pbc'], 
    'CATACAGCGACTACGCATTCTCAT': ['bc26_nbd_mlk', 'bc26_rbk_pbc'], 
    'CGACGGTTAGATTCACCTCTTACA': ['bc27_nbd_mlk', 'bc27_rbk_pbc'], 
    'TGAAACCTAAGAAGGCACCGTATC': ['bc28_nbd_mlk', 'bc28_rbk_pbc'], 
    'CTAGACACCTTGGGTTGACAGACC': ['bc29_nbd_mlk', 'bc29_rbk_pbc'], 
    'TCAGTGAGGATCTACTTCGACCCA': ['bc30_nbd_mlk', 'bc30_rbk_pbc'], 
    'TGCGTACAGCAATCAGTTACATTG': ['bc31_nbd_mlk', 'bc31_rbk_pbc'], 
    'CCAGTAGAAGTCCGACAACGTCAT': ['bc32_nbd_mlk', 'bc32_rbk_pbc'], 
    'CAGACTTGGTACGGTTGGGTAACT': ['bc33_nbd_mlk', 'bc33_rbk_pbc'], 
    'GGACGAAGAACTCAAGTCAAAGGC': ['bc34_nbd_mlk', 'bc34_rbk_pbc'], 
    'CTACTTACGAAGCTGAGGGACTGC': ['bc35_nbd_mlk', 'bc35_rbk_pbc'], 
    'ATGTCCCAGTTAGAGGAGGAAACA': ['bc36_nbd_mlk', 'bc36_rbk_pbc'], 
    'GCTTGCGATTGATGCTTAGTATCA': ['bc37_nbd_mlk', 'bc37_rbk_pbc'], 
    'ACCACAGGAGGACGATACAGAGAA': ['bc38_nbd_mlk', 'bc38_rbk_pbc'], 
    'CCACAGTGTCAACTAGAGCCTCTC': ['bc39_nbd_mlk', 'bc39_rbk_pbc'], 
    'TAGTTTGGATGACCAAGGATAGCC': ['bc40_nbd_mlk', 'bc40_rbk_pbc'], 
    'GGAGTTCGTCCAGAGAAGTACACG': ['bc41_nbd_mlk', 'bc41_rbk_pbc'], 
    'CTACGTGTAAGGCATACCTGCCAG': ['bc42_nbd_mlk', 'bc42_rbk_pbc'], 
    'CTTTCGTTGTTGACTCGACGGTAG': ['bc43_nbd_mlk', 'bc43_rbk_pbc'], 
    'AGTAGAAAGGGTTCCTTCCCACTC': ['bc44_nbd_mlk', 'bc44_rbk_pbc'], 
    'GATCCAACAGAGATGCCTTCAGTG': ['bc45_nbd_mlk', 'bc45_rbk_pbc'], 
    'GCTGTGTTCCACTTCATTCTCCTG': ['bc46_nbd_mlk', 'bc46_rbk_pbc'], 
    'GTGCAACTTTCCCACAGGTAGTTC': ['bc47_nbd_mlk', 'bc47_rbk_pbc'], 
    'CATCTGGAACGTGGTACACCTGTA': ['bc48_nbd_mlk', 'bc48_rbk_pbc'], 
    'ACTGGTGCAGCTTTGAACATCTAG': ['bc49_nbd_mlk', 'bc49_rbk_pbc'], 
    'ATGGACTTTGGTAACTTCCTGCGT': ['bc50_nbd_mlk', 'bc50_rbk_pbc'], 
    'GTTGAATGAGCCTACTGGGTCCTC': ['bc51_nbd_mlk', 'bc51_rbk_pbc'], 
    'TGAGAGACAAGATTGTTCGTGGAC': ['bc52_nbd_mlk', 'bc52_rbk_pbc'], 
    'AGATTCAGACCGTCTCATGCAAAG': ['bc53_nbd_mlk', 'bc53_rbk_pbc'], 
    'CAAGAGCTTTGACTAAGGAGCATG': ['bc54_nbd_mlk', 'bc54_rbk_pbc'], 
    'TGGAAGATGAGACCCTGATCTACG': ['bc55_nbd_mlk', 'bc55_rbk_pbc'], 
    'TCACTACTCAACAGGTGGCATGAA': ['bc56_nbd_mlk', 'bc56_rbk_pbc'], 
    'GCTAGGTCAATCTCCTTCGGAAGT': ['bc57_nbd_mlk', 'bc57_rbk_pbc'], 
    'CAGGTTACTCCTCCGTGAGTCTGA': ['bc58_nbd_mlk', 'bc58_rbk_pbc'], 
    'TCAATCAAGAAGGGAAAGCAAGGT': ['bc59_nbd_mlk', 'bc59_rbk_pbc'], 
    'CATGTTCAACCAAGGCTTCTATGG': ['bc60_nbd_mlk', 'bc60_rbk_pbc'], 
    'AGAGGGTACTATGTGCCTCAGCAC': ['bc61_nbd_mlk', 'bc61_rbk_pbc'], 
    'CACCCACACTTACTTCAGGACGTA': ['bc62_nbd_mlk', 'bc62_rbk_pbc'], 
    'TTCTGAAGTTCCTGGGTCTTGAAC': ['bc63_nbd_mlk', 'bc63_rbk_pbc'], 
    'GACAGACACCGTTCATCGACTTTC': ['bc64_nbd_mlk', 'bc64_rbk_pbc'], 
    'TTCTCAGTCTTCCTCCAGACAAGG': ['bc65_nbd_mlk', 'bc65_rbk_pbc'], 
    'CCGATCCTTGTGGCTTCTAACTTC': ['bc66_nbd_mlk', 'bc66_rbk_pbc'], 
    'GTTTGTCATACTCGTGTGCTCACC': ['bc67_nbd_mlk', 'bc67_rbk_pbc'], 
    'GAATCTAAGCAAACACGAAGGTGG': ['bc68_nbd_mlk', 'bc68_rbk_pbc'], 
    'TACAGTCCGAGCCTCATGTGATCT': ['bc69_nbd_mlk', 'bc69_rbk_pbc'], 
    'ACCGAGATCCTACGAATGGAGTGT': ['bc70_nbd_mlk', 'bc70_rbk_pbc'], 
    'CCTGGGAGCATCAGGTAGTAACAG': ['bc71_nbd_mlk', 'bc71_rbk_pbc'], 
    'TAGCTGACTGTCTTCCATACCGAC': ['bc72_nbd_mlk', 'bc72_rbk_pbc'], 
    'AAGAAACAGGATGACAGAACCCTC': ['bc73_nbd_mlk', 'bc73_rbk_pbc'], 
    'TACAAGCATCCCAACACTTCCACT': ['bc74_nbd_mlk', 'bc74_rbk_pbc'], 
    'GACCATTGTGATGAACCCTGTTGT': ['bc75_nbd_mlk', 'bc75_rbk_pbc'], 
    'ATGCTTGTTACATCAACCCTGGAC': ['bc76_nbd_mlk', 'bc76_rbk_pbc'], 
    'CGACCTGTTTCTCAGGGATACAAC': ['bc77_nbd_mlk', 'bc77_rbk_pbc'], 
    'AACAACCGAACCTTTGAATCAGAA': ['bc78_nbd_mlk', 'bc78_rbk_pbc'], 
    'TCTCGGAGATAGTTCTCACTGCTG': ['bc79_nbd_mlk', 'bc79_rbk_pbc'], 
    'CGGATGAACATAGGATAGCGATTC': ['bc80_nbd_mlk', 'bc80_rbk_pbc'], 
    'CCTCATCTTGTGAAGTTGTTTCGG': ['bc81_nbd_mlk', 'bc81_rbk_pbc'], 
    'ACGGTATGTCGAGTTCCAGGACTA': ['bc82_nbd_mlk', 'bc82_rbk_pbc'], 
    'TGGCTTGATCTAGGTAAGGTCGAA': ['bc83_nbd_mlk', 'bc83_rbk_pbc'], 
    'GTAGTGGACCTAGAACCTGTGCCA': ['bc84_nbd_mlk', 'bc84_rbk_pbc'], 
    'AACGGAGGAGTTAGTTGGATGATC': ['bc85_nbd_mlk', 'bc85_rbk_pbc'], 
    'AGGTGATCCCAACAAGCGTAAGTA': ['bc86_nbd_mlk', 'bc86_rbk_pbc'], 
    'TACATGCTCCTGTTGTTAGGGAGG': ['bc87_nbd_mlk', 'bc87_rbk_pbc'], 
    'TCTTCTACTACCGATCCGAAGCAG': ['bc88_nbd_mlk', 'bc88_rbk_pbc'], 
    'ACAGCATCAATGTTTGGCTAGTTG': ['bc89_nbd_mlk', 'bc89_rbk_pbc'], 
    'GATGTAGAGGGTACGGTTTGAGGC': ['bc90_nbd_mlk', 'bc90_rbk_pbc'], 
    'GGCTCCATAGGAACTCACGCTACT': ['bc91_nbd_mlk', 'bc91_rbk_pbc'], 
    'TTGTGAGTGGAAAGATACAGGACC': ['bc92_nbd_mlk', 'bc92_rbk_pbc'], 
    'AGTTTCCATCACTTCAGACTTGGG': ['bc93_nbd_mlk', 'bc93_rbk_pbc'], 
    'GATTGTCCTCAAACTGCCACCTAC': ['bc94_nbd_mlk', 'bc94_rbk_pbc'], 
    'CCTGTCTGGAAGAAGAATGGACTT': ['bc95_nbd_mlk', 'bc95_rbk_pbc'], 
    'CTGAACGGTCATAGAGTCCACCAT': ['bc96_nbd_mlk', 'bc96_rbk_pbc'], 
    'AAGAAAGTTGTCGGTGTCTTTGTG': ['bc01_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'TCGATTCCGTTTGTAGTCGTCTGT': ['bc02_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'GAGTCTTGTGTCCCAGTTACCAGG': ['bc03_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'TTCGGATTCTATCGTGTTTCCCTA': ['bc04_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'CTTGTCCAGGGTTTGTGTAACCTT': ['bc05_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'TTCTCGCAAAGGCAGAAAGTAGTC': ['bc06_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'GTGTTACCGTGGGAATGAATCCTT': ['bc07_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'TTCAGGGAACAAACCAAGTTACGT': ['bc08_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'AACTAGGCACAGCGAGTCTTGGTT': ['bc09_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'AAGCGTTGAAACCTTTGTCCTCTC': ['bc10_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'GTTTCATCTATCGGAGGGAATGGA': ['bc11_rbk_pbc_pbk_rpb_rab_16s_pcb'], 
    'CAGGTAGAAAGAAGCAGAATCGGA': ['bc12_rbk_pbc_pbk_rab_16s_pcb'], 
    'GTTGAGTTACAAAGCACCGATCAG': ['bc12_rpb']}


#==============================================================================
version = '2025-05-22'  # version of the script
#==============================================================================
def check_version(version):
    try:   
        link = urllib.request.urlopen('https://raw.githubusercontent.com/avierstr/umipore/main/umipore.py').read()
        # find the version-date part of the last version on the webpage
        datepart = re.compile(r'(version.*?)(\d{4}-\d{2}-\d{2})(.*version of the script)')
        x = datepart.search(str(link))
        # the 2nd group of the search is the date
        latest_version = x.group(2)
        # compare the date of this version with the version on the webpage
        if version < latest_version:
            version_name = 'umipore_' + latest_version + '.py' 
            # download latest version
            urllib.request.urlopen('https://raw.githubusercontent.com/avierstr/umipore/main/umipore.py')
            urllib.request.urlretrieve('https://raw.githubusercontent.com/avierstr/umipore/main/umipore.py',
                                        version_name)
            print('\n =====================================================\n'
                  '| NEW VERSION OF umipore AVAILABLE                    |\n'
                  '| https://github.com/avierstr/umipore                 |\n'
                  '| Downloaded latest version as:                       |\n' 
                  '|      ' + version_name + '                          |\n'
                  '| Press ctrl-c to exit                                |\n'
                  ' =====================================================\n')
            t = 10
            while t > 0:
                print('Will continue in ' + str(t) + ' seconds...', end='\r')
                time.sleep(1)
                t -= 1
            # to clear previous line completely   
            print('                                                ', end='\r') 
    except:
        pass
#==============================================================================
def get_arguments():
    
    def range_limited_float_type(arg):
        """ Type function for argparse - a float within some predefined bounds """
        try:
            f = float(arg)
        except ValueError:    
            raise argparse.ArgumentTypeError("Must be a floating point number")
        if f < 0 or f > 1:
            raise argparse.ArgumentTypeError("Argument must be > " + str(0) 
                                             + " and < " + str(1))
        return f

    def valid_file(param): 
        paramlist = [] # make list of files to process
        if os.path.isfile(param): # if input is file
            # check if input file ends on .fastq or .fasta
            base, ext = os.path.splitext(param)
            if ext.lower() not in ('.fasta', '.fastq', '.gz'): 
                raise argparse.ArgumentTypeError('File extension must be .fastq, .fasta or .gz') 
            paramlist.append(param)
        elif os.path.isdir(param): # if input is folder 
            with os.scandir(param) as iterator:
                for file in iterator:
                    if file.name.lower().endswith('.fastq') or file.name.lower().endswith('.fasta')\
                        or file.name.lower().endswith('.gz'):
                        paramlist.append(file.path)
            paramlist.sort()
            if len(paramlist) == 0:
                sys.exit('Can not find files in folder to process.  File extension must be .fastq or .fasta')
        else:
            sys.exit('Can not find a file or folder to process.  File extension must be .fastq or .fasta')
        param = paramlist
        return param

    def dir_path(string):
        string = os.path.join(os.getcwd(), string)
        if not os.path.exists(string):
            os.makedirs(string) # create the folder
        return string
            
    parser = argparse.ArgumentParser(description='Umipore, a tool for demultiplexing nanopore reads',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', required=True, type = valid_file,
                        help='Input folder or file in fastq or fasta format')
    parser.add_argument('-o', '--outputfolder', type=dir_path, 
                         help='Save the results in the specified outputfolder.\n'
                             'Default = folder named as inputfile')
    parser.add_argument('-sf', '--sformat', choices=['auto', 'fasta', 'fastq'], default='auto',
                        help='File format to save the files.  Default = same as inputfile')  
    parser.add_argument('-c', '--compressed', action = 'store_true',
                        help='Compress output files with .gz.')
    parser.add_argument('-min', '--minlength', type = int, default=50,
                        help='Minimum readlenght to process. Default = 50')
    parser.add_argument('-max', '--maxlength', type = int, 
                        help='Maximum readlenght to process.  Default = No limit')
    parser.add_argument('-np', '--nprocesses', type = int, default=4,
                        help='Number of processors to use. Default = 4')  
    parser.add_argument('-sk', '--seq_kit', type = lambda s : s.lower(), 
                        choices = ['lsk', 'lsk109', 'lsk110', 'lsk114', 'cas9', 'rna', 'cdna',
                                   'nbd', 'mlk', 'cs9', 'rad', 'rbk', 'ulk', 'lrk', 'psk', 'pbk',
                                   'rpb', '16s', 'rab', 'pcs', 'dcs', 'pcb'],
                        help='Sequencing kit used in experiment.\n'
                            'Available kits from Nanopore:\n'
                            'LIGATION KITS\n'
                            '  lsk: ligation sequencing kit (lsk109, lsk110, lsk114)\n'
                            '  nbd: native barcoding kit\n'
                            '  mlk: multiplex ligation kit\n'
                            '  cs9: Cas9 sequencing kit             Discontinued\n'
                            'RAPID KITS\n'
                            '  rad: rapid sequencing kit\n'
                            '  rbk: rapid barcoding kit\n'
                            '  ulk: ultra long DNA sequencing kit\n'
                            '  lrk: Field sequencing kit            Discontinued\n'
                            '  psk: PCR sequencing kit              Discontinued\n'
                            '  pbk: PCR barcoding kit               Discontinued\n'
                            '  rpb: rapid PCR barcoding kit\n'
                            '  16s: 16s barcoding kit\n'
                            '  rab: 16s rapid barcoding kit         Discontinued\n'
                            'RNA - cDNA KITS\n'
                            '  rna: direct RNA sequencing kit\n'
                            '  pcs: cDNA PCR sequencing kit\n'
                            '  dcs: direct cDNA sequencing kit      Discontinued\n'
                            '  pcb: cDNA PCR barcoding kit')
    parser.add_argument('-ca', '--custom_adapters', type = str, 
                        help='File with adapters from other sequencing technologies\n'
                            'in case your data is not from Nanopore Sequencing')
    parser.add_argument('-ns', '--no_split', action = 'store_false', default=True,
                        help='Split the reads on middle adapter. Default = True')  

    # barcodes   
    parser.add_argument('-d', '--dual', nargs=2, type = lambda s : s.lower(), 
                        choices = ['nbd', 'pbc', 'custom'], 
                        help='Dual barcodes are used.\n'
                            'Only 2 possible combinations: with the nbd and pbc kit,\n'
                            'or with the nbd kit and custom barcodes'   )  
    # parser.add_argument('-bck', '--bc_kit', type = lambda s : s.lower(), 
    #                     choices = ['nbd', 'mlk', 'pbc', 'rbk',  'pbk', 'rpb', '16s', 'rab', 'pcb'],
    #                     action='append', help='Barcode kit(s) used in experiment.')
    parser.add_argument('-bcn', '--bc_numbers', nargs='+', # action = 'append', 
                        help='Barcode (kit and) numbers used in experiment.\n'
                            'Default = all.\n'
                            'For custom barcodes, it is not using the barcode numbers in\n'
                            'the name but the line numbers in your input file.') 
                        # here possibility to give just numbers, or first kit and than numbers 
                        # nargs='+': multiple args allowed; action='append': multiple times same arg allowed
                        # args are added in nested list
    parser.add_argument('-bcc', '--bc_custom', type = str, 
                        help='File with custom barcodes used in experiment.')
    parser.add_argument('-tr', '--trim', action = 'store_true',
                        help='Trim the barcodes from the reads.  Default = do not trim barcodes')  
    parser.add_argument('-bc_bs', '--bc_both_sides', action = 'store_true', default=True,
                        help='Search for barcodes on both sides.  Default = True')  
    parser.add_argument('-bc_os', '--bc_one_side', action = 'store_true', default=False,
                        help='Search for barcode on one sides.  Default = False')  
    parser.add_argument('-sp', '--search_part', type = int, default=150,
                        help='Part at begin and end of sequence to search for adapters or \n'
                            'barcodes. Default = 150 bp longer than the adapter or barcode')
    parser.add_argument('-er', '--error', type = range_limited_float_type,
                        default=0.15, help='Percentage error allowed in editdistance for\n'
                            'adapters and barcodes. Default = 0.15 (15 percent)')
    # UMIs
    parser.add_argument('-u', '--umi', type = str, 
                        help='File with UMIs used in experiment.')
    parser.add_argument('-umi_bs', '--umi_both_sides', action = 'store_true', default=True,
                        help='Search for UMIs on both sides.  Default = True')  
    parser.add_argument('-umi_os', '--umi_one_side', action = 'store_true', default=False,
                        help='Search for UMIs on one side.  Default = False')  
   

    args = parser.parse_args()
    # print(args)
    # sys.exit()
    return args
#==============================================================================
def save_arguments(): # save all settings in the result.txt file
    outputfolder = args.outputfolder
    with open(os.path.join(outputfolder,'results.txt'), 'a') as rf:
        rf.write('-----------------------------------------------------------\n')
        rf.write('umipore version: ' + version + '\n')  
        rf.write('-----------------------------------------------------------\n')
        rf.write('- date and time = ' + datetime.datetime.now().strftime(
            "%B %d, %Y, %I:%M%p") + '\n')    
        rf.write('- input file = ' + os.path.abspath(infolder_file) + '\n')
        # rf.write('- input file = ' + (','.join(args.input)) + '\n')
        rf.write('- output folder = ' + os.path.normpath(args.outputfolder) + '\n')
        rf.write('- minlength = ' + str(args.minlength) + '\n')
        rf.write('- maxlength = ' + str(args.maxlength) + '\n')
        rf.write('- maxreads = ' + str(args.maxreads) + '\n')
        rf.write('- used_reads = ' + '\n')
        rf.write('- allreads = ' + str(args.allreads) + '\n')
        rf.write('- n_processes = ' + str(args.nprocesses) + '\n')
        rf.write('- similar_genes = ' + str(args.similar_genes) + '\n')
        rf.write('- similar_species_groups = ' + str(args.similar_species_groups) 
                 + '\n')
        rf.write('- similar_species = ' + str(args.similar_species) + '\n')
        rf.write('- similar_consensus = ' + str(args.similar_consensus) + '\n')
        rf.write('- length_difference_consensus = ' + str(args.length_diff_consensus) + '\n')
        rf.write('- save_fastq = ' + str(args.save_fastq) + '\n')
        rf.write('- random = ' + str(args.random) + '\n')
        rf.write('- compare_all = ' + str(args.all) + '\n')
        rf.write('- alignment = ' + str(args.alignment) + '\n')
        rf.write('- ambiguous bases = ' + str(args.ambiguous) + '\n')
        rf.write('- histogram_only = ' + str(args.histogram_only) + '\n')
        rf.write('-----------------------------------------------------------\n')
#==============================================================================
def align(q, t, m, a, k):
    s = edlib.align(q, t, mode=m, task=a, k=k, # task='path' locations
                    additionalEqualities=[("R", "A"), ("R", "G"),
                          ("Y", "C"), ("Y", "T"),
                          ("M", "A"), ("M", "C"),
                          ("K", "G"), ("K", "T"),
                          ("S", "G"), ("S", "C"),
                          ("W", "A"), ("W", "T"),
                          ("N", "A"), ("N", "T"), 
                          ("N", "G"), ("N", "C")])
    return s
#==============================================================================
def compl_reverse(self):
    '''
    make the complement reverse of the sequence
    '''
    inp  = 'ATCGRYKMSWN' # translate table for complement
    outp = 'TAGCYRMKSWN'
    complement = ''.maketrans(inp, outp)
    R = (self[::-1]).translate(complement)  # complement reverse
    return R
#============================================================================== 
def put_in_dict(key, value, dic):
    if key in dic:
        adapt = dic.get(key)
        adapt.append(value)
        dic[key] = adapt
    else:
        dic[key] = [value]
    return dic
#============================================================================== 
def put_in_dict2(key, value, dic):
    if key in dic:
        adapt = dic.get(key)
        adapt.append(value)
        dic[key] = adapt
    else:
        if type(value) is list:
            dic[key] = value
        else:
            dic[key] = [value]
    return dic
#==============================================================================
# ADAPTERS PART
#==============================================================================
def read_adapters(seq_kit):
    '''
    load all adapters, if adapters are the same for multiple kits, only load 
    them once but assign multiple names.
    '''
    # adapters = {}
    # with open(os.path.join(infolder, 'adapters.txt'), 'r') as adap:
    #     for record in SeqIO.parse(adap, "fasta"):
    #         if not record.id.startswith('#'): # comments in adapter file
    #             name = record.id.lower()
    #             seq = str(record.seq).upper() 
    #             if seq in adapters:
    #                 adapt = adapters.get(seq)
    #                 adapt.append(name)
    #                 adapters[seq] = adapt
    #             else:
    #                 adapters[seq] = [name]
    
    # search adapters needed for this analyses
    middle_adap = {}
    front_adap = {}
    rear_adap = {}
    for k in adapters:
        for n in adapters[k]:
            # for m in middle_kit: # find middle adapters
            if n.lower().find('middle') != -1: 
                put_in_dict(k, n, middle_adap)
            if seq_kit is not None:
                if n.lower().find(seq_kit) != -1 and n.lower().find('adapter') != -1 and n.lower().find('front') != -1:
                    put_in_dict(k, n, front_adap)
                if n.lower().find(seq_kit) != -1 and n.lower().find('adapter') != -1 and n.lower().find('rear') != -1:
                    put_in_dict(k, n, rear_adap)
                    
    middle_adap_list = [[middle_adap[k][0], k] for k in middle_adap]
    middle_front_bc = middle_rear_bc = []
    return adapters, front_adap, rear_adap, middle_adap_list, middle_front_bc, middle_rear_bc
#==============================================================================
def read_custom_adap(custom_adap):
    '''
    load custom adapters from file (only needed if other sequencing technologies were
    used like Illumina, Pacbio,...)
    '''
    adapters = {}
    with open(os.path.join(infolder, custom_adap), newline='') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read()) # ckeck if comma or tab separated
        csvfile.seek(0) # got back to begin of file
        reader = csv.reader(csvfile, dialect)
        for row in reader:
            if any(x.strip() for x in row): # remove empty lines
                if not row[0].startswith('#'): # comment in line
                    n = 'adapter_' + row[0].lower() + '_front'
                    k = row[1].strip().replace(' ','').upper()
                    if len(k) > 0:
                        put_in_dict(k, n, adapters) # front adapter
                    n = 'adapter_' + row[0].lower() + '_rear'
                    k = row[2].strip().replace(' ', '').upper()
                    if len(k) > 0:
                        put_in_dict(k, n, adapters) # rear adapter      

    # search adapters needed for this analyses
    middle_adap = {}
    front_adap = {}
    rear_adap = {}
    for k in adapters:
        for n in adapters[k]:
            # for m in middle_kit: # find middle adapters
            if n.lower().find('middle') != -1: 
                put_in_dict(k, n, middle_adap)
            if n.lower().find('adapter') != -1 and n.lower().find('front') != -1 and n.lower().find('middle') == -1:
                put_in_dict(k, n, front_adap)
            if n.lower().find('adapter') != -1 and n.lower().find('rear') != -1 and n.lower().find('middle') == -1:
                put_in_dict(k, n, rear_adap)
                    
    middle_adap_list = [[middle_adap[k][0], k] for k in middle_adap]
    middle_front_bc = middle_rear_bc = []
    return adapters, front_adap, rear_adap, middle_adap_list, middle_front_bc, middle_rear_bc
#==============================================================================
def create_middle_bc_list(used_bc_list):
    '''
    Create a list with degenarated "middle barcodes".
    Convert the list to front and rear middle barcodes
    '''
    middle_dict = {}
    middle_bc_list = []
    for x in used_bc_list:
        put_in_dict(x[1][-16:-1], x[1], middle_dict) # check the last 15 bases
    if len(used_bc_list)*0.2 < len(middle_dict): # if approximate same amount as in list -> no degeneration
        middle_dict = {}
        for x in used_bc_list:
            put_in_dict(x[1][-6:-1], x[1], middle_dict) # check the last 5 bases
    if len(used_bc_list)*0.2 < len(middle_dict): # if approximate same amount as in list -> no degeneration
        middle_dict = {}
        for x in used_bc_list:
            put_in_dict(x[1], x[1], middle_dict) # use them complete
    for i, x in enumerate(middle_dict.values()): 
        # print(x)
        middle_seq = []
        p = -1
        try:
            while True:
                b = set([z[p] for z in x])
                if len(b) > 1:
                    middle_seq.insert(0, 'N')
                else:
                    middle_seq.insert(0, b.pop())
                p -= 1
        except IndexError:
            pass
        middle_bc_list.append([[[' ', 'middle_' + str(i)]], ''.join(middle_seq)])
    middle_front_bc, middle_rear_bc = create_search_list(middle_bc_list, '_')
    return middle_front_bc, middle_rear_bc
#==============================================================================
def filter_locations(locations):
    """
    check if there are multiple locations where the adapter fits at the same
    start position and only use the first for each start position.
    """
    if len(locations) > 1:  
        for l in range(0, len(locations)-1):
            for m in range(l+1, len(locations)):
                if locations[l][0] == locations[m][0]:
                    locations[m] = ('','')
        locations = [x for x in locations if x != ('','')]
    return locations
#============================================================================== 
def get_longest_hit(loc):
    '''
    If searching for dual barcodes with ligation, it is best to search for the single
    barcode without ligation and also for the dual barcode.  If the ligation was not perfect, 
    a pcr fragment with single barcode can be ligated to an other pcr fragment with dual barcode.
    This filters out the longest hit.
    '''
    for l in range(0, len(loc)-1):
        for m in range(l+1, len(loc)):
            for i, k in enumerate(loc[l]):
                for n in loc[m]:
                    z = set(k)
                    z = set([x for x in z if x != ''])
                    z.update(set(n))
                    if len(z) != 4:
                        z = list(z)
                        z.sort()
                        loc[l][i] = ('','')
    loc = [[y] for x in loc for y in x if y != ('','')]
    return loc
#============================================================================== 
def split_middle_keep(record, locations):
    """
    split reads if (a) middle adapter(s) (is) are found and KEEP adapter
    This is important for PCR fragments with barcode or UMI to keep the information in the read
    If compared with rear bc: keep untill end of locations
    If compared with front bc: keep from begin of locations
    
     _!!!!!!!!! You can not include the ligated barcodes to the pcr-primer fragment to search because
    it is possible to ligate 2 pcr fragments together without the ligation barcode !!!
    -> testen of HW mode dit oplost !!
    """
    seqlist = []
    b = 0 # begin

    for x in locations: # parts in the locations
        record2 = record[b:x]
        seqlist.append(record2)
        b = x
    record2 = record[b:] # last part of the sequence
    seqlist.append(record2)
    return seqlist
#============================================================================== 
def process_middle_barcodes(readlist, middle_adapters, middle_front_bc, middle_rear_bc, error, search_part):
    '''
    check if the one of the barcodes/primers are found somewhere in the read.  
    split reads if (a) middle adapter(s) (is) are found and KEEP adapter
    This is important for PCR fragments with barcode or UMI to keep the information in the read
    If compared with rear bc: keep untill end of locations
    If compared with front bc: keep from begin of locations
    
    If min length is given, check if parts meet the min length.
    '''
    min_length = args.minlength
    search_list = [x for x in [middle_adapters, middle_front_bc, middle_rear_bc] if len(x) > 0]
    readlist2 = []
    for record in readlist:
        locations = set()
        seq = [x for x in record.seq] # make a list of the sequence
        seq2 = seq[:]
        d = 0
        for lis in search_list:
            """
            if a 100% hit is found and there is also a 99% hit, it will only show the 100% hit
            for that reason the script has to scan 2 or more times
            """
            k1 = -1 # minimum error
            score = [' ']
            while len(score) > 0:
                score = []
                for n, BC in lis:
                    k = (len(BC)-BC.count('N'))*error # N matches everything, don't count those as length
                    m = 'HW'
                    a = 'locations'
                    s = align(BC, seq2, m, a, k) 
                    if k > (s['editDistance']) > k1: # if a hit is found
                        # print(s['locations'])
                        sl = filter_locations(s['locations'])
                        # print('--' + str(sl))
                        score.append([s['editDistance'], sl])
                if len(score) > 0:
                    if len(score) > 1:
                        score.sort(key=lambda x: x[0])
                        score = [x for x in score if x[0] == score[0][0]] # get the bc with the same editdistance
                    loc = [y for x, y in score] # get all possible locations

                    for j in range(len(loc)):
                        for n in loc[j]:
                            # print(n)
                            begin, end = n
                            seq2[begin:end] = '-'*(end-begin) # remove the hit zone from the intermediate sequence
                    
                    if len(loc) > 1:
                        loc = get_longest_hit(loc)
                    if lis is middle_adapters:
                        loc = [[y[0], y[1]+1] for x in loc for y in x]
                        loc = [y for x in loc for y in x]
                    elif lis is middle_front_bc:
                        loc = [y[0] for x in loc for y in x]
                    elif lis is middle_rear_bc:
                        loc = [y[1]+1 for x in loc for y in x]
                    locations.update(loc)
                k1 += 1
        if len(locations) > 0:
            if len(locations) > 2:
                with middle_split.get_lock():
                    middle_split.value += 1 # count number of seq cut in 3 or more pieces
            locations = list(locations)
            # print(locations)
            locations.sort()
            seqlist = split_middle_keep(record, locations)
            basename = record.id
            # print('---\n>\n' + str(record.seq))
            for record in seqlist:
                if len(record.seq) >= min_length:
                    # print('>\n' + str(record.seq))
                    name = basename + '_' + str(d)
                    record.id = name
                    readlist2.append(record)
                    d += 1
        else:
            readlist2.append(record)
    readlist = readlist2
    return readlist

#============================================================================== 
# def split_middle_remove(record, locations):
#     """
#     split reads if (a) middle adapter(s) (is) are found and REMOVE adapter
#     """
#     seqlist = []
#     b = 0 # begin
#     for x in locations: # parts in the locations
#         record2 = record[b:x[0]]
#         seqlist.append(record2)
#         b = x[1]+1
#     record2 = record[b:] # last part of the sequence
#     seqlist.append(record2)
#     return seqlist
# #============================================================================== 
# def process_middle(readlist, middle_adapters, error):
#     '''
#     check if the middle adapter is found somewhere in the read.  If min length is given,
#     check if parts meet the min length.
#     '''
#     min_length = args.minlength

#     for BC in middle_adapters:
#         readlist2 = []
#         for record in readlist:
#             k = len(BC)*error
#             s = align(BC, record.seq, k) 
#             d = 0
#             if k > (s['editDistance']) > -1: # if a hit is found
#                 locations = filter_locations(s['locations'])
#                 seqlist = split_middle_remove(record, locations)
#                 basename = record.id
#                 for record in seqlist:
#                     with middle_split.get_lock():
#                         middle_split.value += 1 # count number of reads that have minimum one cut
#                     if len(record.seq) >= min_length:
#                         name = basename + '_' + str(d)
#                         record.id = name
#                         readlist2.append(record)
#                         d += 1
#             else:
#                 readlist2.append(record)
#         readlist = readlist2
#     return readlist
#==============================================================================
def remove_front_end_adapters(readlist, front_adap, rear_adap, error, search_part):
    '''
    remove the adapters at the beginning and end of the read if found
    '''
    min_length = args.minlength
    adapter_length = len(list(front_adap.keys())[0]) # get the length of the adapter
    search_part += adapter_length # adjust the search part based on the adapter length
    readlist2 = []
    for record in readlist:
        # fr_trim = 0
        # check for adapter at the beginning of read
        score = []
        fr = str(record.seq)[0:search_part]
        for BC in front_adap:
            k = (len(BC)-BC.count('N'))*error # N matches everythin, don't count those as length
            m = 'HW'
            a = 'locations'
            s = align(BC, fr, m, a, k) 
            if k > (s['editDistance']) > -1: # if a hit is found
                score.append([s['editDistance'], s['locations'][0][1]])
        if len(score) > 0:
            if len(score) > 1:
                score.sort(key=lambda x: x[0])
            fr_trim = score[0][1]+1
            record = record[fr_trim:] # trim front
            with trimmed_front.get_lock():
                trimmed_front.value +=1
        # check for adapter at end of read
        score = []
        er = str(record.seq)[-search_part:]
        for BC in rear_adap:
            k = (len(BC)-BC.count('N'))*error # N matches everythin, don't count those as length
            m = 'HW'
            a = 'locations'
            s = align(BC, er, m, a, k) 
            if k > (s['editDistance']) > -1: # if a hit is found
                score.append([s['editDistance'], s['locations'][0][0]])
        if len(score) > 0:
            if len(score) > 1:
                score.sort(key=lambda x: x[0])
            end_trim = score[0][1]-search_part
            record = record[:end_trim] # trim end
            with trimmed_end.get_lock():
                trimmed_end.value += 1
        if len(record.seq) >= min_length:
            readlist2.append(record)
    return readlist2
#============================================================================== 
# -------------- BARCODES PART ----------------
#==============================================================================
def prep_bc_numbers(bc_numbers, len_bclist):
    '''
    check the given barcode numbers and BC kit, if none given, use all 96; 
    if ranges are given, fill the range.
    '''
    bc_numbers2 = []
    if bc_numbers is not None: # if name or numbers are given
        if type(bc_numbers[0]) is not list: # if it is a nested list, the numbers are already prepared
            sublist = []
            for i, x in enumerate(bc_numbers):
                try: # check if it is a number
                    x = int(x)
                    if sublist[0] == 'custom': # not looking for number but for line number in file
                        '''
                        line    F  -  R
                        line1   1  -  2
                        line2   3  -  4
                        line3   5  -  6
                        line4   7  -  8
                        line5   9  - 10
                        line6   11 - 12
                        '''
                        if x < 0: # if there are negative numbers
                            e = abs(int(bc_numbers[i]))*2  # end number
                            b = int(bc_numbers[i-1])*2 -1 # begin number
                            for y in range(b+1, e+1): # fill range
                                sublist.append(y)
                        else:
                            sublist.append(int(x)*2 -1)  
                            sublist.append(int(x)*2)  
                    else: # normal barcodes: looking for the number
                        if x < 0: # if there are negative numbers
                            e = abs(int(bc_numbers[i])) # end number
                            b = int(bc_numbers[i-1]) # begin number
                            for y in range(b+1, e+1): # fill range
                                sublist.append(y)
                        else:
                            sublist.append(int(x))      
                except ValueError: # it is text
                    if x in ['nbd', 'mlk', 'pbc', 'rbk', 'pbk', 'rpb', '16s', 'rab', 'pcb', 'custom']: # it is the BC kit name
                        if len(sublist) > 0:
                            bc_numbers2.append(sublist)
                            sublist = []
                        sublist.append(x)
                    else:
                        print(str(x) + ' is not a correct value or possible barcode kit.')
                        print('Examples: nbd, mlk, pbc, rbk, pbk, rpb, 16s, rab, pcb, custom')
                        print('Examples: 1 -15 17 18 19 -23')
                        sys.exit()      
            bc_numbers2.append(sublist)
            bc_numbers = bc_numbers2   

        for sl in bc_numbers: # sublist in list
            # print(sl)
            if not type(sl[0]) is str: # if it is not the BC kit name but a number
                if bc_kit in ['nbd', 'mlk', 'pbc', 'rbk', 'pbk', 'rpb', '16s', 'rab', 'pcb']: # it is the BC kit name
                    sl.insert(0, bc_kit) # insert name
                elif bool(bc_custom) is True:
                    sl.insert(0, 'custom') # insert name
                else: # hier nog kijken voor die pbc of custom kit
                    print('Missing a possible barcode kit.')
                    print('Examples: nbd, mlk, pbc, rbk, pbk, rpb, 16s, rab, pcb, custom')
                    sys.exit()    
            if not type(sl[-1]) is int: # if the last item is not a number, use the complete list
                if sl[0] in ['rpb', '16s', 'pcb', 'rab']:
                    sl.extend(x for x in range(1,25))
                if sl[0] in ['custom']:
                    if len_bclist != '_':
                        sl.extend(x for x in range(1, int(len_bclist/2) + 1))
                else:
                    sl.extend(x for x in range(1,97))

    return bc_numbers
#==============================================================================
def read_barcodes(adapters, bc_kit, bc_numbers):
    '''
    Read the barcodes that are used in the experiment.  Depending on the kit, barcodes 
    can be different and can have a other flanking regions (see chemistry technical document 
    on Nanopore website). If 2 kits are used (PBC and NBD) to make dual barcoded
    samples, they are concatenated in this part. 
    Users can enter BC kit name and numbers, or just the numbers if only one kit is used.
    If there is no F or R => use 'f' in "front_bc" and compl-rev as "rear_bc"
    If there is F and R => use 'f' and 'r' in "front_bc" and compl-rev of 'f' and 'r' in "rear_bc"
    
     F-> XXX       crR->yyy
    crF<-xxx            YYY<-R
    
    If there are dual barcodes:
        when only F: connect the F to the previous F and R
        when F and R: connect F to F and R to R (not both ?)
    
      F2->  AAA  F->  XXX            yyy ->crR aaa-> crF2 
    crF2 <- aaa crF<- xxx            YYY <-R   AAA<- F2
    '''     
    bc_numbers = prep_bc_numbers(bc_numbers, '_')
    
    # if bool(dual) is True: # check if dual barcodes are used, for both, bc numbers are needed
    #     if len(bc_numbers) != len(bc_kit) and type(bc_numbers[0]) is list:
    #         print('\nThe barcode numbers are not provided for each barcode kit.')
    #         for x, y in zip_longest(bc_kit, bc_numbers):
    #             print('barcode kit ' + x + ': ' + str(y))
    #         sys.exit('\nPlease give the used barcode numbers for both kits. (-bcn --bc_numbers)')
    #     elif type(bc_numbers[0]) is not list:
    #         print('\nNo barcode numbers for the barcode kits are provided.')
    #         for x in bc_kit:
    #             print('barcode kit ' + x + ': None')
    #         sys.exit('\nPlease give the used barcode numbers for both kits. (-bcn --bc_numbers)')

    # make list with 5'-3' barcodes
    used_bc_list = []

    bc_numbers = [[str(x).zfill(2) for x in sublist] for sublist in bc_numbers] # change 1 to 01, ...
    bc_numbers = sorted(bc_numbers) # make sure the list starts with nbd

    for i, m in enumerate(bc_numbers): # if multiple kits are used sequentially
        r_flank = None
        used_bc_list2 = []
        for ke in adapters:
            for n in adapters[ke]:
                # find flanking regions
                if n.lower().find(m[0]) != -1 and n.lower().find('flank') != -1 \
                    and n.lower().find('forward') != -1 :
                    f_flank = ke
                elif n.lower().find(m[0]) != -1 and n.lower().find('flank') != -1 \
                    and n.lower().find('reverse') != -1 :
                    r_flank = ke
                elif n.lower().find(m[0]) != -1 and n.lower().find('flank') != -1:
                    f_flank = ke # if forw and rev are the same
        # make barcodes with flanking regions
        for num in bc_numbers[i]: # first bc_kit -> first sublist in bc_numbers 
            for ke in adapters:
                for n in adapters[ke]:
                    if n.lower().find(m[0]) != -1 and n.lower().find('bc' + str(num)) != -1: # find barcodes
                        name = 'BC' + num + '_' + m[0].upper()
                        name2 = ['f', name] # put indication 'front' in name 
                        k = f_flank.replace('---BC---', ke)
                        used_bc_list2.append([[name2], k])
                        if r_flank is not None: # if the reverse flank is different from the forward flank 
                            k = r_flank.replace('---BC---', ke)
                            name2 = ['r', name] # put indication 'rear' in name
                            used_bc_list2.append([[name2], k])
        used_bc_list.append(used_bc_list2)
    # print(used_bc_list)
    # sys.exit()   
    # print(used_bc_list)
    used_bc_list = [x for x in used_bc_list if x != []]
    if len(used_bc_list) == 2: # dual barcodes from nanopore 
        # print(used_bc_list)
        used_bc_list = combine_barcodes(used_bc_list[0], used_bc_list[1])
    else:
        used_bc_list = used_bc_list[0] # it is a nested list with one sublist, unnest it
   
    # for x in used_bc_list:
    #     print(x)
    # sys.exit()
    return used_bc_list
#==============================================================================    
def read_custom_barcodes(file, bc_numbers): #hier nog mogelijk maken om bc_numbers te geven
    ''' 
    Read the file with custom barcodes and put them in list.  If "forward" and "reverse" 
    is used in the name, it means they are different, otherwise they are the same.
    '''
    try:
        used_bc_list = []
        with open(os.path.join(infolder, file), newline='') as csvfile:
            dialect = csv.Sniffer().sniff(csvfile.read()) # ckeck if comma or tab separated
            csvfile.seek(0) # got back to begin of file
            reader = csv.reader(csvfile, dialect)
            for row in reader:
                if any(x.strip() for x in row): # remove empty lines
                    if not row[0].startswith('#'): # comment in line
                        name = ['f', row[0]]
                        fbc = row[1].strip().replace(' ','').upper()
                        used_bc_list.append([[name], fbc])
                        name = ['r', row[0]]
                        rbc = row[2].strip().replace(' ', '').upper()
                        used_bc_list.append([[name], rbc])
    except FileNotFoundError:
        print('Can not find custom barcode file ' + file)
        sys.exit()

    bc_numbers = prep_bc_numbers(bc_numbers, len(used_bc_list))  
    
    # with the custom barcodes it is using the line numbers to indicate what barcodes to use
    for sublist in bc_numbers:
        if sublist[0] == 'custom':
            used_bc_list2 = []
            numb = sorted(list(set(sublist[1:]))) # take only the numbers
            for num in numb: # it is the line number
                n = num -1
                used_bc_list2.append(used_bc_list[n])
    used_bc_list = used_bc_list2
    
    return used_bc_list, bc_numbers
#==============================================================================
def combine_barcodes(used_bc_list1, used_bc_list2):
    '''
    If custom (PCR) and kit barcodes (ligation) are combined the sequences of both are
    concatenated here.  "f" and "r" are used to indicate if there are front and rear, or 
    forward and reverse sequences that are different.  If only "f" is used, it means that 
    the rear (reverse) sequence is the same.  "f" and "r" is only important to combine barcodes
    if they have different flanking regions.
    order is 5' - used_bc_list1 - used_bc_list2 - gene - 3'
    '''
    # for x in used_bc_list1:
    #     print(x)

    used1 = set([tuple(x[0] for x in name) for name, seq in used_bc_list1]) # check if f and r are present
    used2 = set([tuple(x[0] for x in name) for name, seq in used_bc_list2]) # check if f and r are present
    
    used_bc_list3 = []
    # print(used_bc_list1)
    if len(used1) == 2 and len(used2) == 2: # if f en r are present in both
        for na, seq in used_bc_list2:
            for name, sequ in used_bc_list1: 
                if na[0][0][0] == 'f' and name[0][0][0] == 'f':
                    name2 = name[:]
                    name2.append(na[0])
                    seq2 = sequ + seq
                    used_bc_list3.append([name2, seq2])
                if na[0][0][0] == 'r' and name[0][0][0] == 'r':
                    name2 = name[:]
                    name2.append(na[0])
                    seq2 = sequ + seq
                    used_bc_list3.append([name2, seq2])
   
    else: # if f en r are not present in both
        for na, seq in used_bc_list2:
            for name, sequ in used_bc_list1: 
               name2 = name[:]
               name2.append(na[0])
               seq2 = sequ + seq
               used_bc_list3.append([name2, seq2])
                
    # elif len(used) == 2 and len(custom) == 1: # if f en r are not present in both
    #     for na, seq in used_custom_bc_list:
    #         for name, sequ in used_bc_list: 
    #             k2 = key + ke
    #             value2 = value[:]
    #             value2.extend(va[:])
    #             used_bc = put_in_dict2(k2, value2, used_bc)
    # else: # only f is present
    #     for na, seq in used_custom_bc_list:
    #         for name, sequ in used_bc_list: 
    #             k2 = key + ke
    #             value2 = value[:]
    #             value2.extend(va[:])
    #             used_bc = put_in_dict2(k2, value2, used_bc)
    # for x in used_bc_list2:
    #     print(x)
    # sys.exit()
    return used_bc_list3
#==============================================================================
def create_search_list(used_bc_list, bc_kit):
    '''
    All bc (also merged ones) are 5'-3' direction.  Clean the name for proper file output
    and make 2 lists used for search: 
    front_bc is for barcodes on the 5' side of each read
    rear_bc is for the complement reverse barcodes on the 3' side of the read.
    Nested lists: [[name, F_bc, R_bc],...]  [[name, F_cr_bc, R_cr_bc],...]
    '''
    bc_kit = ['nbd', 'mlk', 'rbk', 'pbk', 'rpb', '16s', 'rab', 'pcb', 'pbc']
    # if bc_kit is None:
    #     bc_kit = []
    #------------------------
    def remove_kit(name): # remove kit name from total name and remove 'f' and 'r'
        name2 = []
        for _, x in name:
            for y in bc_kit:
                x = x.replace('_' + y.upper(), '')
            name2.append(x)
        name2 = '_'.join(name2)
        return name2
    #------------------------
    used_bc_list = [[remove_kit(name), seq] for name, seq in used_bc_list]
    front_bc = sorted(used_bc_list, key=lambda x: x[0])
    rear_bc = [[name, compl_reverse(seq)] for name, seq in front_bc]
    return front_bc, rear_bc
#==============================================================================
def process_bc_both_sides(readlist, front_bc, error, search_part):
    """
    search for the barcodes on both sides and save the reads in separate files.
    Only search for the 3' barcode that has the same name as on the 5' side. 
    """
    trim = args.trim
    outputfolder = args.outputfolder
    bc_length = len(front_bc[0][1]) # get the length of the barcode
    search_part += bc_length # adjust the search part based on the barcode length
    resultlist = []
    bc_dict = {} # count number of reads per barcode
    MYLOCK = Lock()
    for record in readlist:
        # check for BC in beginning of read
        fr = str(record.seq)[0:search_part]
        scoreF = find_barcode(fr, front_bc, error)
        
        if len(scoreF) > 0: # bc found
            # check for the same BC at the end of the read
            namesF = [x[2] for x in scoreF] # sometimes there are multiple bc with the same distance
            sub_front_bc = [[name, seq] for name, seq in front_bc if name in namesF]# search for same bc in rear_bc
            er = compl_reverse(str(record.seq)[-search_part:])
            scoreR = find_barcode(er, sub_front_bc, error)
            if len(scoreR) > 0:
                namesR = [x[2] for x in scoreR] # sometimes there are multiple bc with the same distance
                z = set(namesF).intersection(set(namesR))
                if len(z) == 1: # same bc found on both sides
                    bc = z.pop()
                else: # len(z) == 0 or len(z) > 1:  f and r are different or multiple possibilities
                    bc = 'Unknown'
            else: # no barcode found on rear side
                bc = 'Unknown'
        else: # no barcode found in front
            bc = 'Unknown'
        if trim == True and bc != 'Unknown':
            fr_trim = [x[1][1] for x in scoreF if x[2] == bc][0] +1
            end_trim = [x[1][1] for x in scoreR if x[2] == bc][0] +1
            record = record[fr_trim:-end_trim] # trim end
        resultlist.append((record, os.path.join(outputfolder, bc)))
        if bc != 'Unknown': # count the number of reads with barcode found
            with barcoded_reads.get_lock():
                barcoded_reads.value += 1
    with MYLOCK:
        for record, name in resultlist:
            save_results(record, name)
            if name in bc_dict:
                bc_dict[name] += 1
            else:
                bc_dict[name] = 1
        resultqueue.put(bc_dict)
#==============================================================================
def process_bc_one_side(readlist, front_bc, error, search_part):
    '''
    search for a barcode (first try 5', if none, try 3') on one side and save the reads 
    in separate files.  If the same bc is found on both sides it is ok, if a different bc
    is found on both sides, put it in Unknown.
    '''
    trim = args.trim
    outputfolder = args.outputfolder
    bc_length = len(front_bc[0][1]) # get the length of the barcode
    search_part += bc_length # adjust the search part based on the barcode length
    resultlist = []
    bc_dict = {} # count number of reads per barcode
    MYLOCK = Lock()
    for record in readlist:
        namesF = []
        namesR = []
        # check for BC in beginning of read
        fr = str(record.seq)[0:search_part]
        scoreF = find_barcode(fr,  front_bc, error)
            
        # check for the same BC at the end of the read
        er = compl_reverse(str(record.seq)[-search_part:])
        scoreR = find_barcode(er, front_bc, error)
        
        if len(scoreF) > 0:
            namesF = [x[2] for x in scoreF] # sometimes there are multiple bc with the same distance
        if len(scoreR) > 0:
            namesR = [x[2] for x in scoreR] # sometimes there are multiple bc with the same distance
        z = set(namesF).intersection(set(namesR))
        if len(z) == 1: # same bc found on both sides
            bc = z.pop()
        elif len(z) > 1: # multiple possibilities
            bc = 'Unknown'
        else: # z == 0
            if len(namesF) == 1 and len(namesR) == 0: # only 1 front bc found
                bc = namesF[0]
            elif len(namesF) == 0 and len(namesR) == 1: # only 1 rear bc found
                bc = namesR[0]
            else: # different front or rear bc found; or no bc
                bc = 'Unknown'
                
        if trim == True and bc != 'Unknown':
            #++++++++++++++++++++++++++
            short_front_bc = [['F', 'GACGACGTTGTAGAGAGTTTGATCMTGGCTCAG'], 
                              ['R', 'GATGGTCGATGACGGTTACCTTGTTACGACTT']]
            #++++++++++++++++++++++++++
            fr_trim = 0
            end_trim = 1 # will be -1 when slicing record below
            if len(scoreF) > 0:
                fr_trim = [x[1][1] for x in scoreF if x[2] == bc][0] +1
            else:
                scoreF = find_barcode(fr,  short_front_bc, error)
                if len(scoreF) > 0:
                    fr_trim = [x[1][1] for x in scoreF][0] +1
            
            if len(scoreR) > 0:
                end_trim = [x[1][1] for x in scoreR if x[2] == bc][0] +1
            else:
                scoreR = find_barcode(er, short_front_bc, error)
                if len(scoreR) > 0:
                    end_trim = [x[1][1] for x in scoreR][0] +1
            record = record[fr_trim:-end_trim] # trim   
            # if len(record2) == 0:
            #     print(str(fr_trim) + '-' + str(end_trim))
            #     print(bc)
            #     print(str(record.seq))
            # print(len(record))
        resultlist.append((record, os.path.join(outputfolder, bc)))
        if bc != 'Unknown': # count the number of reads with barcode found
            with barcoded_reads.get_lock():
                barcoded_reads.value += 1
    with MYLOCK:
        for record, name in resultlist:
            save_results(record, name)
            if name in bc_dict:
                bc_dict[name] += 1
            else:
                bc_dict[name] = 1
        resultqueue.put(bc_dict)
#==============================================================================
def find_barcode(seq, bc_list, error):
    '''
    Find possible barcodes with the same (best) edit distance.  Multiple bc with the
    same edit distance (can also be zero) are stored in the resultlist.
    '''
    score = []
    # tempscore = []
    # fr_trim = 0
    # end_trim = 0
    # bc_name = 'Unknown'
    for item in bc_list:
        name = item[0]
        BC = item[1]
        k = (len(BC)-BC.count('N'))*error # N matches everythin, don't count those as length
        m = 'HW'
        a = 'locations'
        s = align(BC, seq, m, a, k) 
        if k > (s['editDistance']) > -1: # if a approximate hit is found
            # tempscore.append([s, BC, seq, name])
            # print(name)
            # print(s)
            # p = edlib.getNiceAlignment(s, BC, seq)
            # print(p['query_aligned'])
            # print(p['matched_aligned'])
            # print(p['target_aligned'])
            score.append([s['editDistance'], s['locations'][0], name, BC])
    if len(score) > 0:
        if len(score) > 1:
            # for s, BC, seq, name in tempscore:
            #     print(name)
            #     print(s)
            #     p = edlib.getNiceAlignment(s, BC, seq)
            #     print(p['query_aligned'])
            #     print(p['matched_aligned'])
            #     print(p['target_aligned'])
            # print('------')
            score.sort(key=lambda x: x[0])
            score = [x for x in score if x[0] == score[0][0]] # get the bc with the same editdistance
            # print(tempscore)
            # tempscore.sort(key=lambda x: x[0]['editDistance'])
            # print(tempscore)
            # print(score)
            # print('------')
    return score
#==============================================================================
def find_barcode_1(seq, bc_list, error):
    '''
    Find possible barcodes with the same (best) edit distance.  Multiple bc with the
    same edit distance (NOT zero) are stored in the resultlist.  If the edit distance is 
    zero, the search is stopped because this is the best hit.
    '''
    score = []
    for item in bc_list:
        name = item[0]
        BC = item[1]
        k = (len(BC)-BC.count('N'))*error # N matches everythin, don't count those as length
        m = 'HW'
        a = 'locations'
        s = align(BC, seq, m, a, k) 
        if (s['editDistance']) == 0:
            score.append([s['editDistance'], s['locations'][0], name, BC])
            break
        elif k > (s['editDistance']) > -1: # if a approximate hit is found
            score.append([s['editDistance'], s['locations'][0], name, BC])
    if len(score) > 0:
        if len(score) > 1:
            score.sort(key=lambda x: x[0])
            score = [x for x in score if x[0] == score[0][0]] # get the bc with the same editdistance
    return score
#==============================================================================
# UMI PART
#==============================================================================   
def create_alignment(theor_umi, draft_umi):
    # create an alignment out of a list of reads
    alignlist = []
    t = [x for x in theor_umi] # make list of string
    alignlist.append(t)
    number = re.compile(r'\d+') # the numbers to find in edlib result
    symbol = re.compile(r'\D') # the letters or symbols to find in edlib result
    q = [b for b in draft_umi]  
    s = edlib.align(q, t, mode='NW', task='path',
                    additionalEqualities=[("R", "A"), ("R", "G"),
                                          ("Y", "C"), ("Y", "T"),
                                          ("M", "A"), ("M", "C"),
                                          ("K", "G"), ("K", "T"),
                                          ("S", "G"), ("S", "C"),
                                          ("W", "A"), ("W", "T"),
                                          ("N", "A"), ("N", "T"), 
                                          ("N", "G"), ("N", "C")]) # mode must be NW !!
    scorelist = []
    n = re.findall(number, s['cigar'])
    sy = re.findall(symbol, s['cigar']) 
    for x, y in zip(n, sy):
        scorelist += int(x) * [y]
    insertlist = []
    for i, (x, y) in enumerate(zip_longest(t, scorelist)):
        if y == 'I': # insert gap
            insertlist.append(i)
        if y == 'D': # if delete is needed, insert gap in q
            q.insert(i, '-')
    for i in insertlist:
        for z in alignlist: # if insertion in longest sequence, also 
                            # insert in other aligned sequences
            z.insert(i, '-')
    alignlist.append(q)
    return alignlist
#============================================================================== 
'''
1. trim adapters and split middle
2. search all available umis and store reads in temporary separate files (multiprocessing)
3. make search list of umis
4. read the temporary files and filter based on umis
''' 
def read_umis(file):
    ''' 
    Read the file with UMIs and put them in list.  If "forward" and "reverse" 
    is used in the name, it means they are different, otherwise they are the same.
    '''
    try:
        with open(os.path.join(infolder, file), 'r') as inf: # check the fileformat
            line = inf.readline()
            if line[0] == '>':
                fileformat = 'fasta'
            else:
                fileformat = 'csv'
        used_bc_list = []
        if fileformat == 'fasta':
            with open(os.path.join(infolder, file), 'r') as adap:
                for record in SeqIO.parse(adap, "fasta"):
                    if not record.id.startswith('#'): # comments in adapter file
                        name = record.id.lower()
                        bc = str(record.seq).upper() 
                        if name.find('forward') != -1: # forward barcode
                            name = ['f', name.replace('_forward', '').upper()]
                            used_bc_list.append([[name], bc])
                        elif name.find('reverse') != -1: # reverse barcode
                            name = ['r', name.replace('_reverse', '').upper()]
                            used_bc_list.append([[name], bc])
                        else: # if there is no forw or reverse in name
                            name = ['f', name.upper()]    
                            used_bc_list.append([[name], bc])
        else:
            with open(os.path.join(infolder, file), newline='') as csvfile:
                dialect = csv.Sniffer().sniff(csvfile.read()) # ckeck if comma or tab separated
                csvfile.seek(0) # got back to begin of file
                reader = csv.reader(csvfile, dialect)
                for row in reader:
                    if any(x.strip() for x in row): # remove empty lines
                        if not row[0].startswith('#'): # comment in line
                            name = ['f', row[0]]
                            fbc = row[1].strip().replace(' ', '').upper()
                            used_bc_list.append([[name], fbc])
                            name = ['r', row[0]]
                            rbc = row[2].strip().replace(' ', '').upper()
                            used_bc_list.append([[name], rbc])
    except FileNotFoundError:
        print('Can not find UMIs file ' + file)
        sys.exit()
    return used_bc_list
#==============================================================================
def find_umi_sensu_strictu(alignlist):
    ''' 
    Find the theoretical "sensu strictu" umi in the flanking region where Ns are not more than 
    3 bases far from each other (to exclude Ns in the flanks).  Get the part from the 
    alingment and return the real umi 'sensu strictu'.
    AATGATACGGCGACCACCGAGATCNNNYRNNN--YRNNNYRNNNCGACATCGAGGTGCCAAAC
    AATGATACGGCGACCACCGAGATCCGCTAAACGGCACAACGCCTGG-CATTG-GGTGCCAAAC
    
    GTTTGG-C--A-CCT-CG-ATGTCGNNNYRNNNYRNNNYRNNNGATCTCGGTGGTCGCCGTATCATT
    GTTTGGGCGCATCCTTCGGATCT-GGTGCATAACGGTTTGCCCGATCTCGGTGGTCGCCGTATCATT
    '''
    umi = alignlist[0]
    umipos = [] # position of umi in flanking regions
    for i, x in enumerate(umi):
        if x == 'N':
            umipos.append(i)
    a = len(umipos)
    b = len(umipos) - 1
    while b != a : # 
        a = len(umipos)
        for i, x in enumerate(umipos):
            if i == 0: # begin of list
                if umipos[i+1] - umipos[i] > 3:
                    umipos[i] = ''
            if i == len(umipos)-1: # end of list
                if umipos[i] - umipos[i-1] > 3:
                    umipos[i] = ''
        umipos = [x for x in umipos if x != '']
        b = len(umipos)
    begin = umipos[0] 
    end = umipos[-1] + 1
    theo_umi = ''.join(alignlist[0][begin:end]).replace('-', '') # check the theoretical umi without gaps
    exact_umi = ''.join(alignlist[1][begin:end])
    
    if abs(len(theo_umi) - len(exact_umi)) > 1: # if the length difference is more than 1, don't use it
        exact_umi = ''
    return exact_umi
#==============================================================================    
def get_umi_ss(draft_umi, theor_umi):
    '''
    Get the umi 'sensu strictu'' from the sample'
    CTGAGCCAKRATCRAACYCTNNNY-RNNNYRNNNYRNNNATCTCGTATGCCGTCTTCTGCTTG
    CTGAGCCATGATCAAACTCTCCTTTATAGTAATGTAACGATCTCGTATGCCGTCTTCTGCTTG
    '''
    alignlist = create_alignment(theor_umi, draft_umi)
    exact_umi = (find_umi_sensu_strictu(alignlist)).replace('-', '')
    return exact_umi
#==============================================================================
def find_available_umis_os(readlist, front_bc, error, search_part):
    ''' 
    First go over all reads to find all possible front and/or rear umis. 
    Put those in a queue and at the end create a list with all possibilities.  [[front_umi, name], [...]]
    Use those separate names in the umi list.
    Create several lists based on the name.
    Process file by file with the accompanying list.
    '''
    outputfolder = args.outputfolder
    bc_length = len(front_bc[0][1]) # get the length of the barcode
    search_part += bc_length # adjust the search part based on the barcode length
    resultlist = []
    umilist = []
    MYLOCK = Lock()
    for record in readlist:
        namesF = []
        namesR = []
        name = '_'
        # check for UMI in beginning and end of read
        fr = str(record.seq)[0:search_part]
        scoreF = find_barcode(fr, front_bc, error)
        er = compl_reverse(str(record.seq)[-search_part:])
        scoreR = find_barcode(er, front_bc, error)
        
        if len(scoreF) > 0:
            namesF = [x[2] for x in scoreF] # sometimes there are the same UMIs with different names
        if len(scoreR) > 0:
            namesR = [x[2] for x in scoreR] # sometimes there are the same UMIs with different names
        z = set(namesF).intersection(set(namesR))
        if len(z) == 1: # same UMI found on both sides
            name = z.pop()
            with barcoded_reads.get_lock(): # count number of reads with umi
                barcoded_reads.value += 1
        elif len(z) > 1: # multiple possibilities
            name = 'Unknown'
        else: # z == 0
            if len(namesF) == 1 and len(namesR) == 0: # only 1 front bc found
                name = namesF[0]
                with barcoded_reads.get_lock(): # count number of reads with umi
                    barcoded_reads.value += 1
            elif len(namesF) == 0 and len(namesR) == 1: # only 1 rear bc found
                name = namesR[0]
                with barcoded_reads.get_lock(): # count number of reads with umi
                    barcoded_reads.value += 1
            else: # different front and rear umi found; or no UMI
                name = 'Unknown'

        if name != 'Unknown':
            # put front umi in list
            if len(scoreF) > 0:
                scoreF = [x for x in scoreF if x[2] == name]
                begin = scoreF[0][1][0]
                end = scoreF[0][1][1]+1
                draft_umi = fr[begin:end] # found umi with primer flanks
                theor_umi = scoreF[0][3] # theoretical umi with primer flanks
                exact_umi = get_umi_ss(draft_umi, theor_umi) # umi without primer flanks
                if exact_umi != '' :
                    umilist.append([name, exact_umi])
            # put rear umi in list
            if len(scoreR) > 0:
                scoreR = [x for x in scoreR if x[2] == name]
                begin = scoreR[0][1][0]
                end = scoreR[0][1][1]+1
                draft_umi = er[begin:end] # found umi with primer flanks
                theor_umi = scoreR[0][3] # theoretical umi with primer flanks
                exact_umi = get_umi_ss(draft_umi, theor_umi) # umi without primer flanks
                if exact_umi != '' :
                    umilist.append([name, exact_umi])
        filename = 'temp_' + name
        resultlist.append((record, os.path.join(outputfolder, filename)))
    with MYLOCK:
        for record, name in resultlist:
            save_results(record, name)
        umiqueue.put(umilist)
        
        
    # outputfolder = args.outputfolder
    # bc_length = len(front_bc[0][1]) # get the length of the barcode
    # search_part += bc_length # adjust the search part based on the barcode length
    # resultlist = []
    # umilist = []
    # MYLOCK = Lock()
    # for record in readlist:
    #     name = '_'
    #     # check for BC in beginning of read
    #     fr = str(record.seq)[0:search_part]
    #     scoreF = find_barcode(fr, front_bc, error)
    #     if len(scoreF) > 0: # umi found
    #         with barcoded_reads.get_lock(): # count number of reads with umi
    #             barcoded_reads.value += 1
    #         for dis, loc, name in scoreF:
    #             begin = loc[0]
    #             end = loc[1]+1
    #             draft_umi = fr[begin:end]
    #             exact_umi = draft_umi # correct_flanks(draft_umi, name)
    #             if exact_umi != '_':
    #                 umilist.append([name, exact_umi])
    #         # check for the same umi at the end of the read
    #         namesF = [x[2] for x in scoreF] # sometimes there are multiple umis with the same distance
    #         # sub_rear_bc = [[name, seq] for name, seq in rear_bc if name in namesF]# search for same umi in rear_bc
    #     er = str(record.seq)[-search_part:]
    #     scoreR = find_barcode(er, rear_bc, error)
    #     if len(scoreR) > 0:
    #         for dis, loc, name in scoreR:
    #             begin = loc[0]
    #             end = loc[1]+1
    #             draft_umi = compl_reverse(er[begin:end])
    #             exact_umi = draft_umi # correct_flanks(draft_umi, name)
    #             if exact_umi != '_':
    #                 umilist.append([name, exact_umi])
    #         namesR = [x[2] for x in scoreR] # sometimes there are multiple umis with the same distance
    #             # z = set(namesF).intersection(set(namesR))
    #             # if len(z) == 1: # same umi found on both sides
    #                 # name = z.pop()
                   
    #             # else: # len(z) == 0 or len(z) > 1:  f and r are different or multiple possibilities
    #                 # name = 'Unknown'
    #         # else: # no umi found on rear side
    #             # name = 'Unknown'
    #     else: # no barcode found in front
    #         name = 'Unknown'
    #     # if name != 'Unknown':
    #     #     # put front umi in list
    #     #     scoreF = [x for x in scoreF if x[2] == name]
    #     #     begin = scoreF[0][1][0]
    #     #     end = scoreF[0][1][1]+1
    #     #     draft_umi = fr[begin:end]
    #     #     exact_umi = draft_umi # correct_flanks(draft_umi, name)
    #     #     if exact_umi != '_':
    #     #         umilist.append([name, exact_umi])
    #     #     # put rear umi in list
    #     #     scoreR = [x for x in scoreR if x[2] == name]
    #     #     begin = scoreR[0][1][0]
    #     #     end = scoreR[0][1][1]+1
    #     #     draft_umi = compl_reverse(er[begin:end])
    #     #     exact_umi = draft_umi # correct_flanks(draft_umi, name)
    #     #     if exact_umi != '_':
    #     #         umilist.append([name, exact_umi])
    #     filename = 'temp_' + name
       
    #     resultlist.append((record, os.path.join(outputfolder, filename)))
    # with MYLOCK:
    #     for record, name in resultlist:
    #         save_results(record, name)
    #     umiqueue.put(umilist)

#==============================================================================
def find_available_umis_bs(readlist, front_bc, error, search_part):
    ''' 
    First go over all reads to find all possible front umis (find rear umis and make complement-reverse). 
    Put those in a queue and at the end create a list with all possibilities.  [[front_umi, name], [...]]
    If there are barcodes or different umis (16S, operon) involved, first split the reads based on barcodes or umis 
    and save to separate files (dit kan probleem zijn voor 16S en operon !!!!!!).  Use those separate names in the umi list.
    Create several lists based on the name.
    Process file by file with the accompanying list.
   
    Go over the reads again to 
    '''
    outputfolder = args.outputfolder
    bc_length = len(front_bc[0][1]) # get the length of the barcode
    search_part += bc_length # adjust the search part based on the barcode length
    resultlist = []
    umilist = []
    MYLOCK = Lock()
    for record in readlist:
        name = '_'
        # check for BC in beginning of read
        fr = str(record.seq)[0:search_part]
        scoreF = find_barcode(fr, front_bc, error)
        if len(scoreF) > 0: # umi found
            # check for the same umi at the end of the read
            namesF = [x[2] for x in scoreF] # sometimes there are multiple umis with the same distance
            sub_front_bc = [x for x in front_bc if x[0] in namesF]# search for same umi in rear_bc
            er = compl_reverse(str(record.seq)[-search_part:])
            scoreR = find_barcode(er, sub_front_bc, error)
            if len(scoreR) > 0:
                namesR = [x[2] for x in scoreR] # sometimes there are multiple umis with the same distance
                z = set(namesF).intersection(set(namesR))
                if len(z) == 1: # same umi found on both sides
                    name = z.pop()
                    with barcoded_reads.get_lock(): # count number of reads with umi
                        barcoded_reads.value += 1
                else: # len(z) == 0 or len(z) > 1:  f and r are different or multiple possibilities
                    '''Multiple possibilities: sometimes the number of errors is high so that 2 different primers
                    have the same number of errors and none of them are the obvious correct one.  Sometimes the primer
                    on the front and rear side are the same because 2 fragments are ligated together (double in length).'''
                    name = 'Unknown'
            else: # no umi found on rear side
                name = 'Unknown'
        else: # no barcode found in front
            name = 'Unknown'
        if name != 'Unknown':
            # put front umi in list
            scoreF = [x for x in scoreF if x[2] == name]
            begin = scoreF[0][1][0]
            end = scoreF[0][1][1]+1
            draft_umi = fr[begin:end] # found umi with primer flanks
            theor_umi = scoreF[0][3] # theoretical umi with primer flanks
            exact_umi = get_umi_ss(draft_umi, theor_umi) # umi without primer flanks
            if exact_umi != '' :
                umilist.append([name, exact_umi])
            # put rear umi in list
            scoreR = [x for x in scoreR if x[2] == name]
            begin = scoreR[0][1][0]
            end = scoreR[0][1][1]+1
            draft_umi = er[begin:end] # found umi with primer flanks
            theor_umi = scoreR[0][3] # theoretical umi with primer flanks
            exact_umi = get_umi_ss(draft_umi, theor_umi) # umi without primer flanks
            if exact_umi != '' :
                umilist.append([name, exact_umi])
        filename = 'temp_' + name
        resultlist.append((record, os.path.join(outputfolder, filename)))
    with MYLOCK:
        for record, name in resultlist:
            save_results(record, name)
        umiqueue.put(umilist)
#==============================================================================
def process_umi_both_sides(readlist, front_bc, error, search_part):
    ''' 
    Find umis on both sides of the reads.  Save them based on the umi numbers of both sides
    '''
    trim = args.trim
    outputfolder = args.outputfolder
    bc_length = len(front_bc[0][1]) # get the length of the barcode
    search_part += bc_length # adjust the search part based on the barcode length
    resultlist = []
    bc_dict = {} # count number of reads per umi
    MYLOCK = Lock()

    for record in readlist:
        # check for umi in beginning of read
        fr = str(record.seq)[0:search_part]
        scoreF = find_barcode_1(fr, front_bc, error)

        if len(scoreF) > 0: # umi found
            # check for umi at the end of the read
            namesF = [x[2] for x in scoreF] # sometimes there are multiple bc with the same distance
            name = '_'.join(namesF[0].split('_')[:-1]) # base part without number
            numbersF = [x.replace(name + '_', '') for x in namesF] # get the numbers from the front umis
            er = str(record.seq)[-search_part:]
            scoreR = find_barcode_1(er, rear_bc, error)
            if len(scoreR) > 0:
                namesR = [x[2] for x in scoreR] # sometimes there are multiple bc with the same distance
                numbersR = [x.replace(name + '_', '') for x in namesR] # get the numbers from the front umis
                if len(numbersF) > 1 or len(numbersR) > 1:
                    # print(scoreF)
                    # print(scoreR)
                    # print('--')
                    numbersF.sort()
                    numbersR.sort()
                numberlist = sorted([numbersF[0], numbersR[0]]) # sort the numbers
                umi_name = name + '_' + '_'.join(numberlist)
            else: # no barcode found on rear side
                umi_name = 'Unknown'
        else: # no barcode found in front
            umi_name = 'Unknown'
        if trim == True and umi_name != 'Unknown':
            fr_trim = scoreF[0][1][1] +1
            end_trim = scoreR[0][1][0] -search_part
            record = record[fr_trim:end_trim] # trim end
        resultlist.append((record, os.path.join(outputfolder, umi_name)))
        if umi_name != 'Unknown': # count the number of reads with barcode found
            with barcoded_reads.get_lock():
                barcoded_reads.value += 1
    with MYLOCK:
        for record, name in resultlist:
            save_results(record, name)
            if name in bc_dict:
                bc_dict[name] += 1
            else:
                bc_dict[name] = 1
        resultqueue.put(bc_dict)
#==============================================================================
def create_umi_search_list(umiqueue, umi_list, front_bc):
    '''
    Make a list with available umis 'sensu strictu' found in the file.
    '''
    umilist = []
    p = 0
    for umis in iter(umiqueue.get, 'STOP'):# do stuff until infile.get returns "STOP"  
        p += 1
        for name, umi in umis:
            if name != 'Unknown':
                f = 0 # item not found
                for n, l in umilist:
                    if name == n: # if gene is already in list, add UMI to list
                        l.append(umi)
                        f = 1 # item found 
                if f == 0: # item was not found
                    umilist.append([name, [umi]]) # for each gene fragment, create set with UMIs, add to list
     
    for name, ulist in umilist:  # make sublist to reduce processing time
        sub_umi = list(set(ulist))
        pick_name = 'temp_' + name + '.pick'
        with open(os.path.join(outputfolder, pick_name), 'wb') as f:
            pickle.dump(list(set(sub_umi)), f)
#==============================================================================    
def reduce_umis(adapter):
    """
    Compare UMIs with each other to reduce the number of UMIs.  This will reduce the memory
    usage, time needed to find the correct UMI and elimniate sequencing errors.
    """
    umi_todoqueue = Queue(maxsize=1)
    result_umiqueue = Queue()
    progress_count = Value('i', 0) # multiprocessing Value: (interger, 0)
    stopper = Event() # initialize event to signal stop
    nprocesses = args.nprocesses
    
    print('Trying to reduce number of UMIs, this can take some time...')
    d = len(adapter)
    adapter = [[x] for x in adapter]
    # give each UMI a score based on the composition of the bases.  This helps to 
    # reduce the number of comparisons that have to be done.
    for i, x in enumerate(adapter):
        j = 0
        for base in x[0]:
            if base == 'G':
                j += 64
            elif base == 'C':
                j += 16
            elif base == 'T':
                j += 4
            elif base == 'A':
                j += 1
        adapter[i].insert(0, j)
    # sort based on the score
    adapter.sort(key=lambda x: x[0])

    #--------------------------------------------------------------------------  
    def queue_list(adapter): # queue the comparisons in parts for multiprocessing
        l = len(adapter)
        n = 0
        z = 150000 
        k = 0
        while l > z:
            sublist = adapter[n:n+z]
            # umi_todoqueue.put(sublist)
            todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
            with open(todofilename, 'wb', buffering=0) as wf:
                pickle.dump(sublist, wf)
                k += 1
            l = l - z
            n += z
        else:
            sublist = adapter[n:]
            # umi_todoqueue.put(sublist)
            todofilename = os.path.join(outputfolder, 'file_' + str(k) + '.todo')
            with open(todofilename, 'wb', buffering=0) as wf:
                pickle.dump(sublist, wf)
                k += 1
        time.sleep(2)
        
        for dirpath, dirnames, filenames in os.walk(outputfolder):
            filenames = [i for i in filenames if i.endswith('.todo')]
            filenames.sort()
            for name in filenames:
                with open(os.path.join(outputfolder, name), 'rb') as rf:
                    sublist = pickle.load(rf)
                    umi_todoqueue.put(sublist)
                    time.sleep(2)
                    os.remove(os.path.join(outputfolder, name))
        for i in range(nprocesses): # put 'STOP' at the end of the queue for every process
            umi_todoqueue.put("STOP") 
    #--------------------------------------------------------------------------    
    def task():
        process = [Process(target=compare_umis, args=(umi_todoqueue,)) for x in range(nprocesses)]
        for p in process:
            p.start()
        for p in process:
            p.join() 
    #--------------------------------------------------------------------------
    def check_correct_umi(umi):
        A2 = 'NNNYRNNNYRNNNYRNNN'
        k = 0
        m = 'NW'
        a = 'distance'
        s = align(umi, A2, m, a, k)
        dis = (s['editDistance'])   
        return dis
    #--------------------------------------------------------------------------      
    def compare_umis(umi_todoqueue): # compare UMIs with each other base on the score
        for sub_umi in iter(umi_todoqueue.get, 'STOP'): # do stuff until infile.get returns "STOP"
            position = 0
            for position in range(position, len(sub_umi)-1):
                sc = sub_umi[position][0] # score of the UMI
                # if 1 change is allowed for the same length, the score difference can only be the following numbers
                sc_set0 = set([x + sc for x in [3, 12, 15, 48, 60, 63]])
                # if 1 change is allowed for 1 base length difference, the score difference can only be the following numbers
                sc_set1 = set([x + sc for x in [1, 2, 4, 7, 8, 11, 13, 14, 16, 19, 28, 31, 32, 44, 47, 49, 52, 56, 59, 61, 
                                                62, 64, 67, 76, 79, 112, 124, 127]])
                # if 1 change is allowed for 2 bases length difference, the score difference can only be the following numbers
                sc_set2 = set([x + sc for x in [128, 1, 2, 131, 4, 5, 7, 8, 10, 11, 140, 13, 14, 143, 16, 17, 20, 23, 28, 
                                                29, 31, 32, 35, 40, 43, 44, 46, 47, 176, 50, 52, 53, 55, 56, 58, 188, 61, 
                                                62, 191, 65, 68, 71, 77, 80, 83, 92, 95, 113, 116, 125]])
                with progress_count.get_lock():
                    progress_count.value += 1 # count progress
                for position2 in range(position+1, len(sub_umi)):
                    A1 = sub_umi[position][1]
                    if A1 == ' ': 
                        break  # if one of the sequences is removed, no need to compare
                    else:
                        # only compare if the score beween 2 umis is close enough
                        A2 = sub_umi[position2][1]    
                        l = abs(len(A1) - len(A2)) # length difference
                        if l == 0:
                            sc_set = sc_set0
                            msc = sc + 63 # maximum score difference
                        elif l == 1:
                            sc_set = sc_set1
                            msc = sc + 127
                        elif l > 1:
                            sc_set = sc_set2
                            msc = sc + 191
                        else:
                            pass
                        if sub_umi[position2][0] in sc_set: # if score difference is in the set
                            k = 1 + l
                            m = 'NW'
                            a = 'distance'
                            s = align(A1, A2, m, a, k)
                            if k >= (s['editDistance']) > -1: # if a hit is found    
                                if check_correct_umi(A1) == 0: # check if it is a correct umi
                                    sub_umi[position2][1] = ' ' 
                                else:
                                    sub_umi[position][1] = ' ' 
                        elif sub_umi[position2][0] > msc:
                            break # they are sorted, so the next difference will be even bigger      
            sub_umi = [x for x in sub_umi if x[1] != ' ']   
            result_umiqueue.put(sub_umi)
    #--------------------------------------------------------------------------      
    def collect(result_umiqueue, result_adapter): # collect results
        for sub_umi in iter(result_umiqueue.get, 'STOP'): # do stuff until infile.get returns "STOP"
            result_adapter.extend(sub_umi)
    #--------------------------------------------------------------------------      
    def progress(): # show the progress of the camparisons
        print(d)
        while True:
            if stopper.is_set():
                e = '\n'
                f = 100
                print('\r{:.3f} % done'.format(f), end=e)
                break
            else:
                e = ''
                f = round(progress_count.value/(d)*100,3)
            print('\r{:.3f} % done'.format(f), end=e)
            time.sleep(2)
    #--------------------------------------------------------------------------      

    pr = Thread(target=progress)
    pr.start()
    result_adapter = [] # initialize list for thread results
    c = Thread(target=collect, args=(result_umiqueue, result_adapter,))
    c.start()
    # start consumer
    p = Thread(target=task)
    p.start()
    # then start the producer
    q = Thread(target=queue_list, args=(adapter,)) # start reading file(s) in queue
    q.start() 
    # wait for producer to stop
    q.join()
    # wait for consumer to stop
    p.join() # wait until p has finished its work
    result_umiqueue.put('STOP') # put stop in resultqueue when everything is finished
    c.join()
    stopper.set()
    pr.join()
       
    b = len(result_adapter)
    print('Number of UMIs reduced from ' + str(d) + ' to ' + str(b))
    return result_adapter
#==============================================================================
def read_umi_files():
    outputfolder = args.outputfolder
    sformat = args.sformat.lower()
    for dirpath, dirnames, filenames in os.walk(outputfolder):
        for name in sorted(filenames):
            if name.startswith('temp_Unknown'):
                new_name = name.replace('temp_', '')
                os.rename(os.path.join(outputfolder, name), os.path.join(outputfolder, new_name))
            elif name.startswith('temp_') and name.endswith('.pick'):
                # read the umi list corresponding to the file
                nameu = name.replace('temp_', '').replace('.pick', '')
                with open(os.path.join(outputfolder, name), 'rb') as f:
                    umi_list = pickle.load(f)
                    umi_list = reduce_umis(umi_list)
                    # umi_list = [x[1] for x in umi_list]
                    # umi_list = [[[['f', nameu + '_' + str(i)]], umi] for i, [[score], umi] in enumerate(umi_list)]
                    umi_list = [[nameu + '_' + str(i), score, umi] for i, [score, umi] in enumerate(umi_list)]
                    # front_bc, rear_bc = create_search_list(umi_list)
                   
                # os.remove(os.path.join(outputfolder, name))
                # read the umi file corresponding to the umis
                umi_file = name.replace('.pick', '.fastq')
                if os.path.isfile(os.path.join(outputfolder, umi_file)) is False: # check if it is fastq or fasta
                    umi_file = name.replace('.pick', '.fasta')
                if os.path.isfile(os.path.join(outputfolder, umi_file)) is False: # check if it is fastq or fasta
                    break 
                try:
                    with open(os.path.join(outputfolder, umi_file), 'r') as inf: # check the fileformat
                        line = inf.readline()
                        if line[0] == '>':
                            fileformat = 'fasta'
                        elif line[0] == '@':
                            fileformat = 'fastq'
                    if sformat == 'auto':
                        args.sformat = fileformat
                except FileNotFoundError:
                    pass
                stopper.clear() # clear the stop signal for progress
                b = Thread(target=bc_result, args=(resultqueue,))
                b.start()
                p = Thread(target=progress_umi, args=(umi_file, len(umi_list), )) # start progress display
                p.start()
                # first start the consumer
                c = Thread(target=process_umi_queue) # process the reads
                c.start()
                # then start the producer
                q = Thread(target=umi_queue, args=(umi_file, front_bc, rear_bc, )) # start reading file(s) in queue
                q.start() 
                
               
                # wait for producer to stop
                q.join()
                # wait for consumer to stop
                c.join() # wait until c has finished its work
                
                stopper.set() # send stop signal to progress
                p.join() # wait for progress display to stop
                resultqueue.put('STOP') # put stop in resultqueue when everything is finished
                time.sleep(1)
                b.join()

#==============================================================================
def umi_queue(umi_file, front_bc, rear_bc):
    '''
    Make multiprocessing queue to process the umi temporary files
    '''
    fileformat = args.sformat
    readlist = []
    barcoded_reads.value = 0 # initialize values to zero
    file_count.value = 0
    handle = open(os.path.join(outputfolder, umi_file), "r") 
    for record in SeqIO.parse(handle, fileformat):
        record = record.upper() # make all sequences uppercase
        with file_count.get_lock():
            file_count.value += 1 # count total number of reads
        readlist.append(record)
        if file_count.value % 500 == 0:
            todoqueue.put([readlist, front_bc, rear_bc]) 
            readlist = []
            # break
    todoqueue.put([readlist, front_bc, rear_bc]) 
    handle.close()
    os.remove(os.path.join(outputfolder, umi_file))
                
    for i in range(nprocesses): # put 'STOP' at the end of the queue for every process
        todoqueue.put("STOP") 
#==============================================================================
def process_umi_queue():
    '''
    process the reads on multiple cores.  This is a serial process per core.  
    Search for umi's'
    '''
    # outputfolder = args.outputfolder
    # for x in glob.glob(os.path.join(outputfolder, '*')):
    #     os.remove(x)
    nprocesses = args.nprocesses
    search_part = args.search_part
    error = args.error
    def todo(todoqueue):
        for readlist, front_bc, rear_bc in iter(todoqueue.get, 'STOP'):# do stuff until infile.get returns "STOP"
            if bool(umi_bs) is True:
                process_umi_both_sides(readlist, front_bc, error, search_part)
            elif bool(umi_os) is True:
                process_umi_one_side(readlist, front_bc, error, search_part)
                
    try:
        process = [Process(target=todo, args=(todoqueue,)) 
                    for x in range(nprocesses)]
        for p in process:
            p.start()
        for p in process:
            p.join() 
    
    except KeyboardInterrupt:
        print("Shutting processes down")
        # Optionally try to gracefully shut down the worker processes here.
        p.terminate()
        p.join()
        sys.exit()
#==============================================================================
def progress_umi(umi_file, num_umi):  
    '''
    Show the progress of umi precessing.
    Show the result of umi splitting at the end.
    ''' 
    time.sleep(1)
    print('\nProcessing ' + umi_file + ' with ' + str(num_umi) + ' possible UMIs')

    line3 = '|{:^18s}|{:^18s}|'.format('reads', 'reads with')
    line4 = '|{:^18s}|{:^18s}|'.format('read', '2 UMIs found')
    line5 = 39*'-'

    if bc_bs is True:
        print('\nSearching for UMIs on both sides.')
    elif bc_os is True:
        print('\nSearching for UMIs on one side.')
    
    print(line5)
    print(line3)
    print(line4)
    print(line5)

    while True:
        line6 = '\r|{:^18,d}|{:^18,d}|'.format(file_count.value, barcoded_reads.value)
        print(line6, end='')
        if stopper.is_set():
            print(line6, end='\n')
            print(line5)
            print('Reads with UMIs found: {:<16,d}'.format(barcoded_reads.value))
            break
        time.sleep(1)
#==============================================================================
def queue(infile):
    '''
    Make multiprocessing queue to process the reads
    '''
    min_length = args.minlength
    max_length = args.maxlength
    nprocesses = args.nprocesses
    sformat = args.sformat.lower()
    
    # check if it is a gzip file or not, changes the open function
    if infile.endswith('.gz'):
        open_func = gzip.open
    else:  # fasta or fastq
        open_func = open
          
    with open_func(os.path.join(infolder,infile), 'rt') as inf: # check the fileformat
        line = inf.readline()
        if line[0] == '>':
            fileformat = 'fasta'
        elif line[0] == '@':
            fileformat = 'fastq'
    if sformat == 'auto':
        args.sformat = fileformat
    readlist = []
    # file =  'sample_Flo1_BC101.fastq' #"sample_Flo1_BC101.fastq"   #'sample_unk.fastq'
    handle = open_func(os.path.join(infolder,infile), "rt") 
    for record in SeqIO.parse(handle, fileformat):
        record = record.upper() # make all sequences uppercase
        with file_count.get_lock():
            file_count.value += 1 # count total number of reads
        if len(record.seq) >= min_length:
            if max_length is None:
                readlist.append(record)
            else: 
                if len(record.seq) <= max_length:
                    readlist.append(record)
        if file_count.value % 1000 == 0:
            todoqueue.put(readlist) 
            readlist = []
            # break
    todoqueue.put(readlist) 
    for i in range(nprocesses): # put 'STOP' at the end of the queue for every process
        todoqueue.put("STOP")  
    handle.close()
#==============================================================================
def process_queue():
    '''
    process the reads on multiple cores.  This is a serial process per core.  
    First remove adapters if wanted,
    next split reads on middle adapters if wanted,
    next search for barcodes if wanted
    '''
    outputfolder = args.outputfolder
    for x in glob.glob(os.path.join(outputfolder, '*')):
        os.remove(x)
    time.sleep(1)
    middle = args.no_split
    nprocesses = args.nprocesses
    search_part = args.search_part
    error = args.error
    def todo(todoqueue):
        for readlist in iter(todoqueue.get, 'STOP'):# do stuff until infile.get returns "STOP"
            if seq_kit is not None:
                readlist = remove_front_end_adapters(readlist, front_adap, rear_adap, error, search_part)
            if bool(middle) is True:
                readlist = process_middle_barcodes(readlist, middle_adap, middle_front_bc, middle_rear_bc, error, search_part)
                with passed_reads.get_lock():
                    passed_reads.value += len(readlist)# count total number of reads passed length at end
            if bool(bc_bs) is True:
                process_bc_both_sides(readlist, front_bc, error, search_part)
            elif bool(bc_os) is True:
                process_bc_one_side(readlist, front_bc, error, search_part)
            elif bool(umi) is True:
                if bool(umi_bs) is True:
                    find_available_umis_bs(readlist, front_bc, error, search_part)
                elif bool(umi_os) is True:
                    find_available_umis_os(readlist, front_bc, error, search_part)
            else: # save the trimmed reads
                name = os.path.join(outputfolder, infile.replace('.gz', '').replace('.fastq', '_trim').replace('.fasta', '_trim'))
                for record in readlist:
                    save_results(record, name)
                
    try:
        process = [Process(target=todo, args=(todoqueue,)) 
                    for x in range(nprocesses)]
        for p in process:
            p.start()
        for p in process:
            p.join() 
        
    except KeyboardInterrupt:
        print("Shutting processes down")
        # Optionally try to gracefully shut down the worker processes here.
        p.terminate()
        p.join()
        sys.exit()
#==============================================================================
def progress(filename):  
    '''
    Show the progress of trimming the reads from adapters.
    Show the result of barcode splitting at the end.
    '''
    # === adapter information
    middle = args.no_split
    for x in front_adap:
        namelist = []
        for y in front_adap[x]:
            namelist.append(' '.join(y.split('_')[1:2]))
            name = ', '.join(namelist)
        print('Searching for front adapter: ' + str(name).upper() + '\n  (' + str(x) + ') '  )                
    for x in rear_adap:
        namelist = []
        for y in rear_adap[x]:
            namelist.append(' '.join(y.split('_')[1:2]))
            name = ', '.join(namelist)
        print('Searching for rear  adapter: '  + str(name).upper() + '\n  (' + str(x) + ') ' )                  
    if bool(middle) is True:
        search_list = [middle_adap, middle_front_bc, middle_rear_bc]
        for lis in search_list:
            for item in lis:
                name = ''.join(item[0].split('_')[0])
                x = item[1]
                print('Splitting on ' + name + ': '  + '\n  (' + str(x) +')')
            
    # === barcode kit information
    if bool(bc_kit) is True and bool(bc_custom) is False: # if a barcode kit is used
        line1 = "\nSearching for 5' barcodes with flanking regions:"
        line2 = "\nSearching for 3' complement reverse barcodes with flanking regions:"
    # === custom barcode information
    if bool(bc_custom) is True and bool(bc_kit) is False: # if custom barcode is used
        line1 = "\nSearching for 5' custom barcodes:"
        line2 = "\nSearching for 3' complement reverse custom barcodes:"
    # === barcode kit and custom barcode information
    if bool(bc_kit) is True and bool(bc_custom) is True: # if custom and barcode kit are combined
        line1 = "\nSearching for 5' kit and custom barcodes:"
        line2 = "\nSearching for 3' complement reverse kit and custom barcodes:"
        
    # === UMI information
    if bool(umi) is True:
        line1 = "\nSearching for 5' UMIs:"
        line2 = "\nSearching for 3' complement reverse UMIs:"
       
    if bool(bc_kit) is True or bool(bc_custom) is True or bool(umi) is True:
        print(line1)
        for item in front_bc:
            name = item[0]
            x = item[1]
            print(str(name) + '\n  (' + str(x) + ')')
        print(line2)
        for item in rear_bc:
            name = item[0]
            x = item[1]
            print(str(name) + '\n  (' + str(x) + ')')
        # lines to print if there is a search for barcodes or not
        line3 = '|{:^11s}|{:^11s}|{:^11s}|{:^11s}|{:^11s}|'.format('reads', 'trimmed', 'trimmed', 
                                                       'split', 'reads with')
        if bool(umi) is True:
            line4 = '|{:^11s}|{:^11s}|{:^11s}|{:^11s}|{:^11s}|'.format('read', 'front', 'end', 
                                                       'middle', 'UMI')
        else:
            line4 = '|{:^11s}|{:^11s}|{:^11s}|{:^11s}|{:^11s}|'.format('read', 'front', 'end', 
                                                       'middle', 'barcode')
        line5 = 61*'-'
    else:
        line3 = '|{:^11s}|{:^11s}|{:^11s}|{:^11s}|'.format('reads', 'trimmed', 'trimmed', 
                                                       'split')
        line4 = '|{:^11s}|{:^11s}|{:^11s}|{:^11s}|'.format('read', 'front', 'end', 
                                                       'middle')
        line5 = 49*'-'
    
    if bc_bs is True:
        print('\nSearching for barcodes on both sides.')
    elif bc_os is True:
        print('\nSearching for barcodes on one side.')
        
    print('\nProcessing ' + filename)
    print(line5)
    print(line3)
    print(line4)
    print(line5)
    
    while True:
        if bool(bc_kit) is True or bool(bc_custom) is True or bool(umi) is True:
            line6 = '\r|{:^11,d}|{:^11,d}|{:^11,d}|{:^11,d}|{:^11,d}|'.format(file_count.value, trimmed_front.value, 
                                                            trimmed_end.value, middle_split.value, barcoded_reads.value)
        else:
            line6 = '\r|{:^11,d}|{:^11,d}|{:^11,d}|{:^11,d}|'.format(file_count.value, trimmed_front.value, 
                                                            trimmed_end.value, middle_split.value)
        print(line6, end='')
        if stopper.is_set():
            print(line6, end='\n')
            print(line5)
            print('Reads passed minimum length: {:<16,d}'.format(passed_reads.value))
            break
        time.sleep(1)
#==============================================================================
def bc_result(resultqueue):
    '''
    Show the numbers of reads that are found per barcode
    '''
    bc_count = {}
    for bc_dict in iter(resultqueue.get, 'STOP'):# do stuff until infile.get returns "STOP"    
        for k, v in bc_dict.items():
            if k in bc_count:
                bc_count[k] += v
            else:
                bc_count[k] = v
    print('\nNumber of reads per barcode found:')   
    for k, v in sorted(bc_count.items(), key=lambda item: item[1], reverse=True):
        infolder, infile = os.path.split(os.path.realpath(k))
        print('{:30s} {:.>16,d}'.format(infile, v))
#==============================================================================    
def save_results(record, name):
    '''
    save the results in separate files
    '''
    compressed = args.compressed
    if compressed is True: # output files compressed or not
        open_func = gzip.open
    else:  # fasta or fastq
        open_func = open
        
    sformat = args.sformat.lower()
    if sformat == 'fastq':
        name = name + '.fastq'
    else:
        name = name + '.fasta'
    if compressed is True:
        name = name + '.gz'
    with open_func(name, 'at') as writer: 
        try:
            SeqIO.write(record, writer, sformat)
        except ValueError: # can not save fastq from fasta file
            SeqIO.write(record, writer, 'fasta')
#============================================================================== 
if __name__ == '__main__':
    try:
        args = get_arguments()
        check_version(version)
        outputfolder = args.outputfolder
        infolder_file_list = args.input
        nprocesses = args.nprocesses
        seq_kit = args.seq_kit # sequencing kit used
        custom_adap = args.custom_adapters
        init_sformat = args.sformat.lower() # format to save file
        bc_numbers = args.bc_numbers
        stopper = Event() # initialize event to signal stop
        todoqueue = Queue(maxsize = 2)#nprocesses + 1) # max number in queue
        resultqueue = Queue() # queue to store number of reads in BC_files
        umiqueue = Queue() # queue to store the possible umis
        bc_custom = args.bc_custom
        dual = args.dual
        
        # catch incomplete info from users
        if seq_kit == 'nbd' and bool(bc_custom) is True: # check if seq kit is a BC kit and custom is used
            if dual is None:
                print('ERROR:')
                print('Seems you are using dual barcodes: "nbd and custom"')
                print('Please use the "-d or --dual" option for that.')
                sys.exit()
                
        elif seq_kit in ['nbd', 'mlk', 'rbk', 'pbk', 'rpb', '16s', 'rab', 'pcb']: # check if seq kit is a BC kit 
            bc_kit = seq_kit
        
        if bool(dual) is True: # dual barcodes are used
            bc_kit = seq_kit
            if all(x in ['nbd', 'pbc'] for x in dual): # check if both are present
                bc_kit = dual
            if 'custom' in dual and bool(bc_custom) is False:
                print('ERROR:')
                print('You did not mention the custom BC file with the "-bcc" option')
                sys.exit()
            if bc_numbers is not None: # if name of numbers are given
                if not any(x in ['custom', 'nbd', 'pbc'] for x in bc_numbers): # check if one is mentioned
                    print('ERROR:')
                    print('It is not clear for which kit the barcode numbers are used.')
                    sys.exit()
                if 'nbd' not in bc_numbers:
                    bc_numbers.append('nbd')
                if 'custom' in bc_numbers and bool(bc_custom) is False:
                    print('ERROR:')
                    print('You did not mention the custom BC file with the "-bcc" option')
                    sys.exit()
            else:
                bc_numbers = dual

        # if seq_kit in ['nbd', 'mlk', 'rbk', 'pbk', 'rpb', '16s', 'rab', 'pcb']: # check if seq kit is a BC kit 
        #     bc_kit = seq_kit.split()
        #     if bc_numbers is None: # if no name or numbers are given
        #         bc_numbers = []
        #         bc_numbers.append(seq_kit) # use the name of the BC kit in a list
        #         if bool(bc_custom) is True:
        #             bc_numbers.append('custom')
        #     else: # name or numbers are given
        #         if bool(bc_custom) is True: # 2 BC kits used
        #             if not 'custom' and 'nbd' in bc_numbers:
        #                 print('')
        
        # wanneer is bc_numbers nodig ?  als er 2 bc kits samen gebruikt zijn met :
        #     nbd en pbc
        #     nbd en custom
        # anders is er telkens maar 1 mogelijkheid van naam en aantal BC
        
           
        # else: # the seq kit is not a BC kit
        #     if bc_numbers is not None and bool(bc_custom) is False: # if name or numbers are given but no custom BC
        #         bc_kit = ['pbc'] # PCR barcoding expansion kit is the only non-sequencing BC kit
        #     if bool(bc_custom) is True:
        #         bc_numbers = ['custom']
        #         bc_kit = False # no nanopore BC kit used
        #     else:
        #         bc_kit = False
        # print(bc_kit)
        # print(bool(bc_kit))
        
        bc_os = args.bc_one_side
        bc_bs = args.bc_both_sides
        umi = args.umi
        umi_os = args.umi_one_side
        umi_bs = args.umi_both_sides

        if bool(bc_os) is True:
            bc_bs = False
        if bool(bc_kit) is False and bool(bc_custom) is False: # if no barcodes are given
            bc_bs = False
            bc_os = False
        """Dit zal problemen geven als combinatie umi en barcode mogelijk is.  ofwel bc_bc en bc_os site ook gebruiken voor umis ?
        """
             
        for infolder_file in infolder_file_list:
            infolder, infile = os.path.split(os.path.realpath(infolder_file))
            outfolder, ext = os.path.splitext(infile.replace('.gz', ''))
            if not outputfolder: # if outputfolder is not given
                try:
                    args.outputfolder = os.path.join(infolder, outfolder)
                    os.makedirs(args.outputfolder)
                except FileExistsError:
                    pass
        adapters, front_adap, rear_adap, middle_adap, middle_front_bc, middle_rear_bc = read_adapters(seq_kit)
        if bool(custom_adap)is True:
           adapters, front_adap, rear_adap, middle_adap, middle_front_bc, middle_rear_bc = read_custom_adap(custom_adap)
        # --------- BARCODE PART -----------------
        if bool(bc_kit) is True and bool(bc_custom) is False: # if a barcode kit is given
            used_bc = read_barcodes(adapters, bc_kit, bc_numbers)
        if bool(bc_custom) is True and bool(bc_kit) is False: # if a custom barcode file is given
            used_bc, bc_numbers = read_custom_barcodes(bc_custom, bc_numbers)
        if bool(bc_kit) is True and bool(bc_custom) is True: # if kit and custom is given
            used_custom_bc, bc_numbers = read_custom_barcodes(bc_custom, bc_numbers)
            used_bc = read_barcodes(adapters, bc_kit, bc_numbers)
            used_bc = combine_barcodes(used_bc, used_custom_bc)
            
        # -----------UMI PART --------------------
        if bool(umi) is True:
            used_bc = read_umis(umi) # read the file with theoretical UMI's 
        """"
        hier nog combine barcodes bijvoegen: umi en barcode aan geligeerd
        """
        if bool(bc_kit) is True or bool(bc_custom) is True or bool(umi) is True: # if barcodes are used, make list to search front and rear of a read
            middle_front_bc, middle_rear_bc = create_middle_bc_list(used_bc) # if BC are used, create list with middle adapters
            front_bc, rear_bc = create_search_list(used_bc, bc_kit)
        

        # if bool(umi) is True:
        #     find_umi_sensu_strictu(front_bc, rear_bc)
        #     flank_list = find_umi_flanks(front_bc) # make a dict of the flanking sequences of the umis
        
        # Initialize values
        args.sformat = init_sformat
        file_count = Value('i', 0) # multiprocessing Value: (interger, 0) read from file
        trimmed_front = Value('i', 0) # multiprocessing Value: (interger, 0) reads trimmed front
        trimmed_end = Value('i', 0) # multiprocessing Value: (interger, 0) reads trimmed back
        middle_split = Value('i', 0) # multiprocessing Value: (interger, 0) reads split middle
        barcoded_reads = Value('i', 0) # multiprocessing Value: (interger, 0) read with barcode
        passed_reads = Value('i', 0) # multiprocessing Value: (interger, 0) read passed length criteria
        
        # save_arguments() # write all settings in the results.txt file
        stopper.clear() # clear the stop signal for progress
        p = Thread(target=progress, args=(infile,)) # start progress display
        p.start()
        time.sleep(1)
        
        if bool(bc_bs) is True or bool(bc_os) is True:
            b = Thread(target=bc_result, args=(resultqueue,))
            b.start()
        if bool(umi) is True:
            # make sublists of UMIs for each sequenced gene
            umi_list = []
            u = Thread(target=create_umi_search_list, args=(umiqueue, umi_list, front_bc, )) #dit geeft memory problemen
            u.start()
            
        # start the consumer
        c = Thread(target=process_queue) # process the reads
        c.start()
        # start the producer
        q = Thread(target=queue, args=(infile,)) # start reading file(s) in queue
        q.start() 

        q.join() # wait for producer to finish
        c.join() # wait for consumer to finish

        stopper.set() # send stop signal to progress
        p.join() # wait for progress display to stop
        if bool(bc_bs) is True or bool(bc_os) is True:
            resultqueue.put('STOP') # put stop in resultqueue when everything is finished
            b.join()
        
        if bool(umi) is True:
            # bc_bs = True
            umiqueue.put('STOP') # put stop in umiqueue when everything is finished
            u.join()
            
            read_umi_files()

            
    except KeyboardInterrupt:
        sys.exit()

"""
te doen:
- adapters in de file zelf zetten, in goed leesbare lijst -> OK
- barcodes vergelijken tussen de verschillende kits en ook in een goede lijst zetten -> OK in lijst van adapters
- testen om zelf middel barcodes te zoeken uit lijst > OK
- is de middle adapter ook als F en R te vinden ? -> enkel als F -> OK
- chemistry guide nakijken, -> OK
- zorgen dat andere adapters ook kunnen geladen worden -> OK
- automatisch laten zoeken naar adapters en barcodes met random functie om ook reads op einde file te lezen -> niet zo interessant
  door overlap van BC en adapters tussen verschillende kits en dual barcodes
- NBD barcode flanking region nakijken -> OK (is identiek aan Complement reverse)
- alle adapters en benamingen van kits nakijken en controleren -> OK
- argparse deftig uitlijnen -> OK
- mogelijkheid om gz files te lezen en schrijven -> OK
- testen of split on middle ook de primers aan de rand afknipt
- moeten andere adapters ook als middle beschouwd worden ? cas9: idem als lsk; RNA: nakijken bij info kits front en rear adapters 
cDNA PCR sequencing kit   (contains UMI ???)
- nog RNA data zoeken om te testen
- save arguments aanpassen

-nbd: native barcoding kit: als je die selecteerd, moet je geen bc kit meer opgeven, het moet automatisch weten welke bc 
type gebruikt is.

dual barcoding is mogelijk met de pbc gevolgd door nbd (lsk ??) kit.  of met eigen pcr en dan nbd.

UMI:
- gz support voorzien
gecheckt voor nbd en pbc -> ok

python3 umipore_2025-06-04.py -i Unknown.fasta -sk nbd -d custom nbd -bcc custom.csv -bcn custom 1 4 -6 nbd 3 -6
python3 umipore_2025-06-04.py -i Unknown.fasta -sk nbd -d pbc nbd -bcn pbc 1 -6 nbd 3 -6

testen voor alle andere kits

"""