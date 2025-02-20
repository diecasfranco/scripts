#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 14:41:39 2020

@author: Andy.Vierstraete@ugent.be

Download a database from NCBI
Reduce Nucleotide or protein database from NCBI for faster local Blast searches 
by doing size selection and/or removing duplicate sequences
"""

from Bio import SeqIO
import sys
import socket
import time
import datetime
import os
import pickle
from ftplib import FTP
import shutil
import argparse
import tarfile
import hashlib
from multiprocessing import Process, Lock, Queue
from threading import Thread
import edlib
import glob
import subprocess
import re
import urllib.request

genes = [['COI', 'CO1', 'COX1', 'COXI', 'cytochrome oxidase I',  'cytochrome oxidase 1', 
         'mitochondrion complete genome'],
         ['COII', 'CO2', 'COX2', 'COXII', 'cytochrome oxidase 2', 'cytochrome oxidase II', 
         'mitochondrion complete genome'],
         ['COIII', 'CO3', 'COX3', 'COXIII','cytochrome oxidase III', 'cytochrome oxidase 3',
         'mitochondrion complete genome'],
         # nog nakijken voor regex:
         ['12S', 'small subunit ribosomal'],
         ['16S', '16S rRNA', '16S ribosomal RNA',  'small subunit ribosomal', '16S rDNA'],
         ['18S', '18S rRNA', '18S ribosomal', 'small subunit ribosomal', 'SSU rRNA'],
         ['23S', 'large subunit ribosomal'],
         ['26S', 'large subunit ribosomal'],
         ['28S', 'LSU', 'large subunit ribosomal', 'large subunit rRNA','large subunit of rRNA'], 
         ['ITS', 'internal transcribed spacer'],
         ['NADH1', 'nad1', 'NADH dehydrogenase subunit 1', 'nd1',
                  'NADH dehydrogenase subunit I', 'NDI', 'nicotinamide adenine dinucleotide dehydrogenase subunit 1'],
         ['NADH2','nad2', 'NADH dehydrogenase subunit 2', 'ND2', 'NADH-ubiquinone oxidoreductase chain 2'],
         ['psbA', 'photosystem II'],
         ['psaA', 'photosystem I'],
         ['rbcL', 'ribulose 1,5 bisphosphate carboxylase/oxygenase large'],
         ['tufA', 'elongation factor Tu'],
         ['rpb2', 'RNA polymerase II', 'RNA polymerase second', 'rpbII'],
         ['ef1', 'TEF', 'elongation factor 1', 'elongation factor alpha', 'EF-1']
         ]
unknown = ['unidentified', 'unknown', 'unspecified', 'untyped', 'ungrouped', 
           'undetermined', 'undescribed', 'uncultured', 'uncultivated', 'unclassified']
#     soil, sludge, sediment, salal, rumen, \
# root, rhizosphere, rape rhizosphere, rainbow trout, psychrophilic, \
# prokaryote enrichment, primary endosymbiont, nitrogen fixing, \
# mycorrhizal, mucus, mixed, methylotrophic, methanotrophic, methanotroph, \
# methanogenic, metal-contaminated, mesophilic, mercury-resistant, marine, \
# maize, low G+C, leaf litter, iron-reducing, intestinal, ice core, \
# humic, halophilic, haloarchaeon, haloalkaliphilic, groundwater, \
# grassland, glacial, gamma proteobacterium, fungal, fuel, freshwater, \
# fossil, forest, foliar, filamentous, fern, extreme, eukaryote, \
# eubacterial, ericoid, epsilon, epacrid, environmental, enrichment, \
# endosymbiont, endophytic, endocytic, ectomycorrhizal, earthworm, \
# drinking, denitrifying, delta, cyanobacterium enrichment, \
# crenarchaeote enrichment, candidate division, beta proteobacterium, \
# barley, bacterium enrichment, archaeon enrichment, anammox bacterium enrichment, \
#     alpha proteobacterium enrichment, activated sludge, actinobacterium enrichment)]

version = '2024-07-08'  # version of the script
#==============================================================================
def check_version(version):
    try:   
        link = urllib.request.urlopen('https://github.com/avierstr/reduceblastdb/blob/main/reduceblastdb.py').read()
        # find the version-date part of the last version on the webpage
        datepart = re.compile(r'(version.*?)(\d{4}-\d{2}-\d{2})(.*version of the script)')
        x = datepart.search(str(link))
        # the 2nd group of the search is the date
        latest_version = x.group(2)
        # compare the date of this version with the version on the webpage
        if version < latest_version:
            version_name = 'reduceblastdb_' + latest_version + '.py' 
            # download latest version
            urllib.request.urlopen('https://raw.githubusercontent.com/avierstr/reduceblastdb/main/reduceblastdb.py')
            urllib.request.urlretrieve('https://raw.githubusercontent.com/avierstr/reduceblastdb/main/reduceblastdb.py',
                                        version_name)
            print('\n =====================================================\n'
                  '| NEW VERSION OF reduceblastdb AVAILABLE              |\n'
                  '| https://github.com/avierstr/reduceblastdb           |\n'
                  '| Downloaded latest version as:                       |\n' 
                  '|      ' + version_name + '                    |\n'
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
#------------------------------------------------------------------------------
def arguments():
    parser = argparse.ArgumentParser(description='Reduceblastdb: script to reduce NCBI \
    databases to specific taxIDs, length and optional only unique sequences')
    
    subparsers = parser.add_subparsers(dest='subparser_name', help='There are 2 sub-commands \
                                       "tax" (taxonomy) and "db" (database). \
                            Taxonomy is to get information about the NCBI classification. \
                            Database is to select parts of the blast databases to create a \
                                smaller local database.')
    
    taxparser = subparsers.add_parser('tax', help='type "tax --help" for more information about \
                                      taxonomy options')
    taxparser.add_argument('-of', '--outfolder', required=False, 
                        help='Name of the output folder to save the taxonmomy files.  \
                            If none is given, current working directory is used.')
    taxparser.add_argument('-tl', '--taxidlist', required=False, nargs='*',
                        help='Create a TaxID list from NCBI with name, taxid and rank \
                            in alfabatic order (fast). List is in csv format. \
                            By default, everything below genus level is excluded \
                            to limit the file size.  You can enter a taxonomic level to \
                            change the default or use "all" to use all levels. \
                            [kingdom, phylum, class, order, family, genus, species, all] \
                            You can limit the list by giving a TaxID or scientific name')
    taxparser.add_argument('-tt', '--taxtree', required=False, nargs='*',
                        help='Create a "Taxonomic tree" from NCBI with rank, name and taxid \
                            in alfabatic order per taxonomic level (very slow).  Information is \
                            saved in csv file.  This is not a fylogenetic tree file. \
                            By default, everything below genus level is excluded \
                            to limit the file size.  You can enter a taxonomic level to \
                            change the default or use "all" to use all levels. \
                            [kingdom, phylum, class, order, family, genus, species, all] \
                            You can limit the tree by giving a TaxID or scientific name')
    taxparser.add_argument('-l', '--lineage', required=False,
                        help='Create the full taxonomic lineage from a TaxID or scientific name')
    
    dbparser = subparsers.add_parser('blastdb', help='type "db --help" for more information about \
                                      blastdatabase options')
    dbparser.add_argument('-db', '--database', required=True,
                        help='Which database ?  nt, nr, LSU_eukaryote_rRNA, ...')
    dbparser.add_argument('-dl', '--download', required=False, action = 'store_true', 
                        help='Download the needed database from NCBI')
    dbparser.add_argument('-o', '--outfile', required=False, 
                        help='Name of the outputfile saved in the outputfolder "reduced"')
    dbparser.add_argument('-of', '--outfolder', required=False, 
                        help='Name of the output folder to save the database and create custom \
                            databases.  Custom created databases will be saved in the subfolder \
                            "reduced".  If none is given, current working directory is used.')
    dbparser.add_argument('-min', '--minlength', type = int, required=False,
                        help='Minimum sequence length')
    dbparser.add_argument('-max', '--maxlength', type = int, required=False,
                        help='Maximum sequence length')
    dbparser.add_argument('-s', '--select', required=False, nargs='+',
                        help='Select TaxIDs to make a sub-database, use "1" for all species \
                            use one or more TaxIDs or scientific names for one or more groups')
    dbparser.add_argument('-u', '--unknown', required=False, action = 'store_true', 
                        help='Exclude non identified entries (unidentified, unknown, ...)')                       
    dbparser.add_argument('-g', '--gene', required=False, nargs='+',
                        help='Select genes to make a sub-database, \
                            use one or more gene names that are used in the "sequence title" in GenBank') 
    dbparser.add_argument('-a', '--append', required=False, action = 'store_true',
                        help='Append sequences to a previous created fasta file. \
                            (it adds the sequences to the "selected.fasta" file)')                         
    dbparser.add_argument('-r', '--reduce', required=False, action = 'store_true', 
                        help='Reduce the number of sequences by removing highly similar sequences \
                            from each species')
    dbparser.add_argument('-np', '--nproc', type = int, default= 4, required=False,
                        help='Number of processors to use in the "reduce" option.  Default= 4')
    dbparser.add_argument('-mbd', '--makeblastdb', required=False, action = 'store_true', 
                        help='Create the new Blast database from the selected or reduced file')
    # parser.add_argument('-id', '--identity', type = int,  required=False,
    #                     help='Identity between sequences.  Default= 97')
    args = parser.parse_args()
    return args
#------------------------------------------------------------------------------
def save_arguments(): # save all settings in the result.txt file
    outpfile = args.outfile
    try:
        with open(os.path.join(outfolder, outpfile + '_results.txt'), 'r') as rf:
            l = rf.readlines()
        with open(os.path.join(outfolder, outpfile + '_results.txt'), 'w') as rf:
            rf.writelines(l)
    except FileNotFoundError:  
        with open(os.path.join(outfolder, outpfile + '_results.txt'), 'w') as rf:
            rf.write('-----------------------------------------------------------\n')
            rf.write('reduceblastdb version: ' + version + '\n')  
            rf.write('-----------------------------------------------------------\n')
            rf.write('- date and time = ' + datetime.datetime.now().strftime("%B %d, %Y, %I:%M%p") + '\n')  
            rf.write('- database = ' + args.database + '\n')
            rf.write('- download = ' + str(args.download) + '\n')
            rf.write('- output file = ' + args.outfile + '\n')
            rf.write('- minlength = ' + str(args.minlength) + '\n')
            rf.write('- maxlength = ' + str(args.maxlength) + '\n')
            rf.write('- select TaxIDs = ' + str(args.select) + '\n')
            rf.write('- reduce = ' + str(args.reduce) + '\n')
            rf.write('- n_processes = ' + str(args.nproc) + '\n')
            rf.write('- makeblastdb = ' + str(args.makeblastdb) + '\n')
#------------------------------------------------------------------------------        
def diskspace(fsize):  #check if there is enough diskspace left
    total, used, free = shutil.disk_usage(dlfolder)
    freeGB = round(free/(1024*1024*1024))
    fsizeGB = round(fsize/(1024*1024*1024))
    print(str(fsizeGB) + ' GB to download, ' + str(freeGB) + ' GB free on disk')
    if freeGB < fsizeGB*3:
        print('!!! Only ' + str(freeGB) + ' GB free on disk, approximately ' + str(fsizeGB*3) +
              ' GB is needed !!!')
        sys.exit('not enough disk space free')
#==============================================================================  
def distance(X1,X2):  # calculate the similarity of 2 sequences
    '''
    this is to compare 2 sequences where the shortest is more than 92% similar with the longest
    in its whole length. (updated 2024_06_23: contained error)
    '''
    iden = 0
    if len(X1) > len(X2): # check which one is longer
        A2 = X1
        A1 = X2
    else:
        A1 = X1
        A2 = X2
    if len(A1) > 100000:  # split long sequences in parts to save time and memory
        part = 500 # fast check to see if there is a similar part
        A1x = A1[0:part]
        k = len(A1x)*0.08 # max 8% difference
        s = edlib.align(A1x, A2, task='distance', mode='HW', k=k) 
        distance = s['editDistance']
        if distance == -1: # not similar
            iden = 0
        else: 
            part = 30000 # is fastest (tested between 5.000 en 30.000)
            b1 = b2 = 0 # begin for A1 and A2
            distotal = 0
            while b1 <= len(A1):
                A1x = A1[b1:b1+part]
                A2x = A2[b2:]
                k = len(A1x)*0.08 # max 8% difference
                s = edlib.align(A1x, A2x, task='locations', mode='HW', k=k) 
                distance = s['editDistance']            
                if distance == -1: # not similar
                    distotal = len(A1)
                    break
                else:
                    distotal += distance
                    b1 += part
                    try:
                        b2 = s['locations'][0][1] # begin for A2
                    except IndexError:
                        b2 = 0
            iden = round(1 - distotal/len(A1),7)  
    else:
        k = len(A1)*0.08 # max 8% difference
        s = edlib.align(A1, A2, task='distance', mode='HW', k=k)
        distance = s['editDistance']
        if distance >= 0: 
            iden = round(1 - distance/len(A1),3)
    return iden
#==============================================================================
def tax_rank(items):
    check_for_files() 
    print('Building database')
    taxdict = {}
    with open(os.path.join(dlfolder, 'names.dmp'), 'r') as f: # read file with taxid and scientific names
        for line in f:
            x = line.strip('\t|\n').split('\t|\t')
            if x[3] == 'scientific name':
                taxdict[x[0]] = x[1]  # taxID : scientific name
                
    with open(os.path.join(dlfolder, 'nodes.dmp'), 'r') as n: # read file with taxID, parent node and rank
        for line in n:
            x = line.strip('\t|\n').split('\t|\t')
            sci_name = taxdict.get(x[0])
            taxdict[x[0]] = [sci_name, x[1], x[2]] # taxID : [scientific name, parent, rank]
    del taxdict['1'] # remove the root
    
    rank = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'all']
    rankorder = ['superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'superclass',
            'class', 'subclass', 'infraclass', 'cohort', 'superorder', 'order', 'suborder',
            'superfamily', 'family', 'subfamily', 'genus', 'species']

    delset = {'genotype', 'serogroup', 'pathogroup', 'varietas', 'morph', 'species subgroup',  
        'subspecies', 'forma specialis', 'forma', 'isolate', 'biotype', 'serotype', 'strain'}

    # check if a rank is given as parameter
    d = list(set(items).intersection(set(rank)))
    if len(d) == 1:
        d = d[0]
        if d == 'all':
            print('Including all taxonomic levels')
            delset = {}
        else: 
            print('Excluding everything below ' + d + ' level')
            i = rankorder.index(d)
            delset.update(rankorder[i+1:])
    elif len(d) > 1:
        print('Only 1 taxonomic level can be used')
        sys.exit()
    else: # nothing is given
        d = 'genus' # default value 
        print('Excluding everything below ' + d + ' level')
        i = rankorder.index(d)
        delset.update(rankorder[i+1:])
        
    # check if a taxonomic group is given
    t = list(set(items).difference(set(rank)))
    if len(t) == 1:
        t = t[0]
        if str(t).isalpha(): # Check if input is name
            for key, value in taxdict.items():
                if value[0].lower() == str(t).lower():
                    t = key
        try:  # catch errors in search name or number
            _ = taxdict[t][0]
        except KeyError:
            if str(t).isalpha():
                print('"' + str(t) + '"' + ' is not an NCBI official scientific name')
            else: 
                print('"' + str(t) + '"' + ' is not an NCBI official taxID')
            sys.exit()
        print('Searching for ' + str(taxdict[t][2]) + ' ' + str(taxdict[t][0]) + ' with taxid ' + str(t))
    elif len(t) > 1:
        print('Only 1 taxID or scientific name can be used')
        sys.exit()
    else: # nothing is given
        t = '1'
        print('Starting at the root')
    return delset, taxdict, t
#==============================================================================
def create_taxidlist(items): 
    delset, taxdict, t = tax_rank(items) # deleteset, taxdict, taxid

    taxlist = taxdict.items()
    # remove taxons from delset
    taxlist = [[key, value] for key, value in taxlist if value[2] not in delset]
    # build a parent dict to increase search speed
    parentdict = {}
    for key, value in taxlist:
        parent = value[1]
        if value[1] in parentdict: # parent
            kids = parentdict.get(parent)
            kids.add(key) # add the kid
            parentdict[parent] = kids
        else:
            parentdict[parent] = {key} # create set
            
    taxlist = []
    tempparent = []
    parent = [t]
    name, _, rank = taxdict.get(t)
    taxlist.append([t, name, rank])
    m = 0
    # start at parents, look for kids
    while len(parent) != 0:  
        print('Searching level ' + str(m))
        for p in parent:
            sublist = []
            kids = parentdict.get(p)
            if kids is not None:
                for k in kids:
                    name, _, rank = taxdict.get(k)
                    # taxlist.append([k, name, rank])
                    sublist.append([k, name, rank])
                # sublist.sort(key=lambda x: x[2].lower()) #sort list based on name
                # try:
                #     i = [x[1] for x in treelist].index(p) + 1
                #     treelist[i:i] = sublist # insert sublist in list with list slicing
                # except ValueError:
            taxlist.extend(sublist)
            tempparent.extend([x[0] for x in sublist])
        parent = tempparent[:]
        tempparent = []
        m += 1
    
    taxlist = sorted(taxlist, key=lambda x: x[1])
    with open(os.path.join(outfolder,'taxidlist.csv'), 'w') as taxidfile:
        for k, name, rank in taxlist:
            print(name + ',' + str(k) + ',' + rank, file=taxidfile)
    print('Scientific names, TaxID and rank saved in file ' + os.path.join(outfolder,"taxidlist.csv"))
#==============================================================================
def create_taxtree(items): 
    delset, taxdict, t = tax_rank(items) # deleteset, taxdict, taxid

    # build a parent dict to increase search speed
    parentdict = {}
    taxlist = taxdict.items()
    # remove taxons from delset
    taxlist = [[key, value] for key, value in taxlist if value[2] not in delset]
    for key, value in taxlist:
        parent = value[1]
        if value[1] in parentdict: # parent
            kids = parentdict.get(parent)
            kids.add(key) # add the kid
            parentdict[parent] = kids
        else:
            parentdict[parent] = {key} # create set
       
    treelist = []
    tempparent = []
    parent = [t]
    m = 0
    name, _, rank = taxdict.get(t)
    treelist.append([m, t, name, rank])
    # start at parents, look for kids
    while len(parent) != 0:  
        print('Searching level ' + str(m))
        for p in parent:
            sublist = []
            kids = parentdict.get(p)
            if kids is not None:
                for k in kids:
                    name, _, rank = taxdict.get(k)
                    sublist.append([m, k, name, rank])
                sublist.sort(key=lambda x: x[2].lower()) #sort list based on name
                try:
                    i = [x[1] for x in treelist].index(p) + 1
                    treelist[i:i] = sublist # insert sublist in list with list slicing
                except ValueError:
                    treelist.extend(sublist)
                tempparent.extend([x[1] for x in sublist])
        parent = tempparent[:]
        tempparent = []
        m += 1

    with open(os.path.join(outfolder,'taxidtree.csv'), 'w') as taxidfile:
        for x in treelist:
            if x[3] == 'kingdom':
                i = '\u25BA'  # hexa symbol
            elif x[3] == 'phylum':
                i = '\u25CF'
            elif x[3] == 'class':
                i = '\u25CB'
            elif x[3] == 'order':
                i = '\u2022'
            elif x[3] == 'family':
                i = '\u25AA'
            elif x[3] == 'genus':
                i = '*'
            elif x[3] == 'species':
                i = '+'
            else:
                i = '-'
            print('   '*x[0] + i + ' ' + x[3].title() + ',' + x[2] + ',' + x[1], file=taxidfile)
    print('Rank, Scientific names and TaxID saved in file ' + os.path.join(outfolder,"taxidtree.csv"))
#==============================================================================
def create_lineage(taxid):
    items = [] # items for function below has to be a list
    items.append(taxid)
    delset, taxdict, t = tax_rank(items) # deleteset, taxdict, taxid

    # search for parents in the nodes
    taxid = t
    parents = []
    while taxid != '1':
        parents.append(taxid)
        taxid = taxdict[taxid][1]
    print('{:15s} {:30s} {:>20s}\n'.format('Rank', 'Scientific name', 'TaxID'))
    for x in parents[::-1]:
        print('{:15s} {:30s} {:>20s}'.format(taxdict[x][2].title(), taxdict[x][0], x))
#==============================================================================
# 1. DOWNLOAD BLAST DATABASE FILES FROM INTERNET
# download ftp://ftp.ncbi.nlm.nih.gov/blast/db/
#==============================================================================
def get_files(filename):
    def count_files(filename):
        try:
            ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=3600)
            ftp.login("anonymous", "")
            ftp.cwd("/blast/db/")
            x = list(ftp.mlsd(facts=['size', 'modify'])) # list the files
            x = [y for y in x if y[0].split('.')[0] == filename] # select only wanted files
            fsize = sum([int(y[1]['size']) for y in x]) # check filesize needed to download
            diskspace(fsize) # check is there is enough space on disk
            db_version = datetime.datetime.strptime(x[0][1]['modify'], '%Y%m%d%H%M%S')
            filelist = [y[0] for y in x if y[0].endswith('.gz')]
            try:
                with open(os.path.join(dlfolder,'database_version.txt'), 'r') as f: # check database version
                    s_version = f.readline().strip()
                    if str(s_version) == str(db_version): # same version, maybe files have already been downloaded
                        downloaded_files = f.readlines() # read which files have been downloaded
                        # print(downloaded_files)
                        downloaded_files = [x.strip() for x in downloaded_files]
                        # print(downloaded_files)
                        filelist_update = sorted(list(set(filelist).difference(set(downloaded_files))))
                        filelist = filelist_update
                    else:
                        try:  # remove previous file if exists
                            for x in glob.glob(os.path.join(dlfolder, filename + '*')):
                                # if not x.endswith('.tar.gz'):
                                os.remove(x)
                            time.sleep(1)
                        except FileNotFoundError:
                            pass
                        with open(os.path.join(dlfolder,'database_version.txt'), 'w') as f: # store database version
                            print(db_version, file=f)
            except FileNotFoundError:
                try:  # remove previous file if exists
                    for x in glob.glob(os.path.join(dlfolder, filename + '*')):
                        # if not x.endswith('.tar.gz'):
                        os.remove(x)
                    time.sleep(1)
                except FileNotFoundError:
                    pass
                with open(os.path.join(dlfolder, 'database_version.txt'), 'w') as f: # store database version
                    print(db_version, file=f)
            ftp.quit()
        except socket.timeout:
            count_files(filename) # restart after timeout
        return filelist
    #--------------------------------------------------------------------------
    def download(ftp, file):
        filesize = ftp.size(file) # check size of file to be downloaded
        lfile = os.path.join(dlfolder, file) # local file
        try:
            endpoint = os.path.getsize(lfile) # check size if part is already downloaded
        except FileNotFoundError:
            endpoint = 0
        while endpoint < filesize: # restart downloading as long as it is not complete
            if os.path.exists(lfile):  # check if a part is already downloaded
                endpoint = os.path.getsize(lfile)  # where did download stop ?
                print('Downloading rest of ' + file)
                # resume downloading where it was stopped
                ftp.retrbinary('RETR ' + file , open(lfile, 'ab').write, rest=endpoint, blocksize=33554432)
                endpoint = os.path.getsize(lfile)  # where did download stop ?
            else: # if it has not been downloaded yet
                print('Downloading ' + file)
                ftp.retrbinary('RETR ' + file , open(lfile, 'wb').write, blocksize=33554432)
                endpoint = os.path.getsize(lfile)  # where did download stop ?
    #--------------------------------------------------------------------------
    def verify(downloadedqueue):
        for filen in iter(downloadedqueue.get, 'STOP'): # do stuff until infile.get returns "STOP"
            lfilen = os.path.join(dlfolder, filen) # local file
            print('verifying ' + filen)
            with open(lfilen, 'rb') as f:
                file_hash = hashlib.md5()
                chunk = f.read(128*1000)
                file_hash.update(chunk)
                while chunk:
                    chunk = f.read(128*1000)
                    file_hash.update(chunk)
            with open(lfilen + '.md5', 'r') as f:
                if (f.read().split()[0]) == file_hash.hexdigest():
                    print('--> ' + filen + ' verified and OK !')
                    print('--> Decompressing ' + filen )
                    my_tar = tarfile.open(lfilen, mode='r:gz')  # 'r|gz' data stream of blocks
                    if filen == 'taxdump.tar.gz':
                        my_tar.extract('names.dmp', dlfolder)
                        my_tar.extract('nodes.dmp', dlfolder)
                    else:
                        my_tar.extractall(dlfolder)
                    my_tar.close()
                    os.remove(lfilen)
                    os.remove(lfilen + '.md5')
                    with open(os.path.join(dlfolder, 'database_version.txt'), 'a') as f: # store downloaded filenames
                        print(filen, file=f)
                    print('--> ' + filen + ' decompressed')
                else:
                    print(filen + ' is corrupt.  Removing file and download again.')
                    os.remove(lfilen)
                    os.remove(lfilen + '.md5')
                    with open(os.path.join(dlfolder,'corrupt.txt'), 'a') as f: # store corrupted files in file for re-download
                        print(filen, file=f)
    #--------------------------------------------------------------------------
    def ftp_process(filelist):
        try:
            for file in filelist:
                ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=3600) # increase buffersize to 32 MB for ncbi  33554432
                ftp.login("anonymous", "")
                ftp.cwd("/blast/db/")
                filemd5 = file + '.md5'
                download(ftp, filemd5)
                download(ftp, file)
                downloadedqueue.put(file) # make queue of downladed files for decompression
                ftp.quit()   
        except socket.timeout:
            ftp_process() # restart after timeout
    #--------------------------------------------------------------------------
    def ftp_taxonomy(file):
        try:
            ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=3600) # increase buffersize to 32 MB for ncbi  33554432
            ftp.login("anonymous", "")
            ftp.cwd("/pub/taxonomy/")
            filemd5 = file + '.md5'
            download(ftp, filemd5)
            download(ftp, file)
            downloadedqueue.put(file) # make queue of downladed files for decompression
            ftp.quit()   
        except socket.timeout:
            ftp_process() # restart after timeout
    #--------------------------------------------------------------------------
    def download_verify_decompress(filelist):
        v = Thread(target=verify, args=(downloadedqueue,)) # start verify/decompress
        v.start()
        ftp_taxonomy('taxdump.tar.gz') # donwload taxonmy files, needed later
        d = Thread(target=ftp_process, args=(filelist,)) # start downloading
        d.start() 
        d.join() # wait for downloading to finish
        downloadedqueue.put('STOP') # add stop to verify/decompress queue
        v.join() # wait for verify/decompress to finish
    #--------------------------------------------------------------------------
    try:    
        os.remove(os.path.join(dlfolder, 'corrupt.txt')) 
    except FileNotFoundError:
        pass
    
    downloadedqueue = Queue()
    filelist = count_files(filename)
    download_verify_decompress(filelist)
    
    try:
        with open(os.path.join(dlfolder, 'corrupt.txt'), 'r') as rf: # if corrupted files, try again
            filelist = rf.readlines()
            os.remove(os.path.join(dlfolder, 'corrupt.txt'))
            download_verify_decompress(filelist)
    except FileNotFoundError:
        pass
#==============================================================================
# 2. SELECT THE TAXIDs YOU WANT TO USE
# You can make a selection on the taxid to reduce the database size
#==============================================================================   
def select_sequences(taxids, outpfile):
    check_for_files()
    print('Loading taxID information from names.dmp and nodes.dmp')
    taxdict = {}
    with open(os.path.join(dlfolder,'names.dmp'), 'r') as f: # read file with taxid and scientific names
        for line in f:
            x = line.strip('\t|\n').split('\t|\t')
            if x[3] == 'scientific name':
                taxdict[x[0]] = x[1]  # taxID - scientific name
    with open(os.path.join(dlfolder,'nodes.dmp'), 'r') as n: # read file with taxid, parent node and rank
        for line in n:
            x = line.strip('\t|\n').split('\t|\t')
            sci_name = taxdict.get(x[0])
            taxdict[x[0]] = [sci_name, x[1], x[2]] # taxID - [scientific name, parent, rank]
    
    taxidlist = []
    for y in taxids: # list of taxids
        if str(y).isalpha(): # Check if input is name
            for key, value in taxdict.items():
                if value[0].lower() == str(y).lower():
                    y = key
            taxidlist.append(y)    
        else: # input is a number
            taxidlist.append(y)
    
    print('Searching all species for:')
    with open(os.path.join(outfolder, outpfile + '_results.txt'), 'a') as rf:
        print('Searching all species for:', file=rf)
    for tax in taxidlist:
        try:  # catch errors in search name or number
            print('- ' + str(taxdict[tax][2]) + ' ' + str(taxdict[tax][0]) + 
            ' with taxid ' + str(tax))
            with open(os.path.join(outfolder, outpfile + '_results.txt'), 'a') as rf:
                print('- ' + str(taxdict[tax][2]) + ' ' + str(taxdict[tax][0]) + 
                ' with taxid ' + str(tax), file=rf)
        except KeyError:
            if str(tax).isalpha():
                print('"' + str(tax) + '"' + ' is not an NCBI official scientific name')
                sys.exit()
            else: 
                print('"' + str(tax) + '"' + ' is not an NCBI official taxID')
                sys.exit()
        
    # search for childeren in the nodes
    species = set()
    for tax in taxidlist:
        childs = set()
        for key, value in taxdict.items():
            if value[1] == str(tax): # initial position in search
                childs.add(key)
    
        b =  len(childs)
        e = b + 1
        while e > b:  # keep searching for descendants until no more found
            b = e
            for key, value in taxdict.items():
                if value[1] in childs:
                    childs.add(key)
            e = len(childs)
    
        subspecies = set()
        for key in taxdict: # select all species (nodes, genera, family are not sequences)
            if key in childs:
                if 'species' in taxdict[key][2]:
                    subspecies.add(key)
                    # print(taxid[key])
        print(str(len(subspecies)) + ' species found for ' + str(taxdict[tax][0]))
        with open(os.path.join(outfolder, outpfile + '_results.txt'), 'a') as rf:
            print(str(len(subspecies)) + ' species found for ' + str(taxdict[tax][0]), file=rf)
        species.update(subspecies)

    specieslist = sorted(list(species))
    with open(os.path.join(outfolder, 'taxidlist.txt'), 'w') as tl:
        for x in specieslist:
            print(x, file=tl)
    
    # Extract accession numbers, taxid and lengths for the sequences 
    print('Searching all accessions, taxids and lengths from the ' + filename + ' database for:')
    for tax in taxidlist:
        print('- ' + str(taxdict[tax][2]) + ' ' + str(taxdict[tax][0]))
       
    command = ['blastdbcmd']
    command.append('-db')
    command.append(os.path.join(dlfolder, filename))
    command.append('-taxidlist')
    command.append(os.path.join(outfolder, 'taxidlist.txt'))
    # command.append('-target_only')
    command.append('-out')
    command.append(os.path.join(outfolder, outpfile + '_selected_taxid.txt'))
    command.append('-outfmt')
    command.append('%a %T %l %i %t')
    # print(command)
    result = subprocess.run(command, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0 :
        print('----> Error ' + str(result.stderr) + '\n')
    else:
        print('----> Selection finished !')  
    
    size_select(outpfile) # size select if wanted, or prepare seqid file
                
    print('Extracting all sequences from:') 
    for tax in taxidlist:   
        print('- ' + str(taxdict[tax][2]) + ' ' + str(taxdict[tax][0]) + 
          ' with taxid ' + str(tax) + ' from the ' + filename + ' database')
    
    command = ['blastdbcmd']
    command.append('-db')
    command.append(os.path.join(dlfolder, filename))
    command.append('-entry_batch')
    command.append(os.path.join(outfolder, outpfile + '_seqid.txt'))
    command.append('-out')
    command.append(os.path.join(outfolder, outpfile + '_selected.fasta'))
    # command.append('-target_only')
    # print(command)
    result = subprocess.run(command, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0 :
        print('----> Error ' + str(result.stderr) + '\n')
    else:
        print('----> Extraction finished !')  
    
    remove_duplicates(outpfile) # blastdbcmd extracts duplicate sequences that are identical 
    
    os.remove(os.path.join(outfolder, 'taxidlist.txt'))
    os.remove(os.path.join(outfolder, outpfile + '_seqid.txt'))
#==============================================================================
# 3. REMOVE DUPLICATES AND DO SIZE SELECTION IF WANTED
# If user only want sequences of a certain size to reduce the database
#==============================================================================
def check_db_type(outpfile):
    # check the type of the database: nucleotide or protein
    # prot = ['A',      'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T',      'V', 'W', 'Y']
    # nucl = ['A', 'B', 'C', 'D',           'G', 'H',      'K',      'M', 'N',           'R', 'S', 'T', 'U', 'V', 'W', 'Y']
    dbset = set() # to check if nucleotide or protein sequence
    t = 0
    inputfile = open(os.path.join(outfolder, outpfile + '_selected.fasta'), "r") 
    for record in SeqIO.parse(inputfile, "fasta"):
        t += 1
        # check if nucleotide or protein
        if t < 10:
            seq = record.seq
            if len([letter for letter in seq.upper() if letter in "EFILPQ"]) > 0: # is protein
                dbset.add('prot')
            else:
                dbset.add('nucl')
    inputfile.close()
    if len(dbset) > 1:
        print('Not clear if database is nucleotide or protein')
        sys.exit()
    else: 
        dbtype = list(dbset)[0]
        with open(os.path.join(outfolder, outpfile + '_dbtype.tmp'), 'w') as dbt:
            print(dbtype, file= dbt)
    return dbtype
#--------------------------------------------------------------------------
def remove_duplicates(outpfile): # blastdbcmd extracts duplicate sequences that are identical
    print('Searching for duplicates based on accession number')
    duplicateset = set() 
    dupfile = os.path.join(outfolder, outpfile + '_no_duplicates.fasta')
    try:
        os.remove(dupfile)    
    except FileNotFoundError:
        pass
    try:
        inputfile = open(os.path.join(outfolder, outpfile + '_selected.fasta'), "r") 
    except FileNotFoundError:
        print('File not found')
    
    _ = check_db_type(outpfile)
    
    t = 0
    for record in SeqIO.parse(inputfile, "fasta"):
        t += 1
        if record.id not in duplicateset: # check for duplicates 
            duplicateset.add(record.id)
            with open(dupfile, 'a') as outputfile:
                SeqIO.write(record, outputfile, "fasta")
    inputfile.close()
    print(str(t - len(duplicateset)) + ' duplicate sequences removed')
    print(str(len(duplicateset)) + ' unique sequences selected ')
    with open(os.path.join(outfolder, outpfile + '_results.txt'), 'a') as rf:
        print(str(t - len(duplicateset)) + ' duplicate sequences removed', file=rf)
    os.remove(os.path.join(outfolder, outpfile + '_selected.fasta'))
    os.rename(dupfile, 
              os.path.join(outfolder, outpfile + '_selected.fasta'))
#--------------------------------------------------------------------------            
def size_select(outpfile):
    seqidfile = os.path.join(outfolder, outpfile + '_seqid.txt')
    gene = args.gene
    t = 0 # total
    s = 0 # size select
    g = 0 # gene
    u = 0 # unknown
    if gene is not None: # search for seq with specific genes
        geneset = set()
        for y in genes: # list of variations in gene names
            for x in gene: # list of wanted genes
                if any(x.upper() == z.upper() for z in y):
                    geneset.update(y) # list of gene names to search for
                else:
                    geneset.add(x)
        print('Selecting sequences with the following terms in the title:')
        for x in sorted(list(geneset)):
            print('  ' + x)
    
    if args.unknown is True: 
        print('Excluding sequences when the title starts with:')
        for x in sorted(unknown):
            print('  ' + x)
        
    if args.minlength is not None or args.maxlength is not None: # do size selection
        print('Selecting sequences between the size limits...')
        outfile = os.path.join(outfolder, outpfile + '_sizeselected.txt')
        minlen = args.minlength
        maxlen = args.maxlength
        with open(os.path.join(outfolder, outpfile + '_selected_taxid.txt'), 'r') as ss:
            for line in ss:
                line = line.strip()
                _, _, l, seqid, title = line.split(' ', 4)
                t += 1 #counting number of sequences scanned
                
                # check if it is the correct gene
                if gene is not None:
                    if any(re.search(x.replace(' ', '.*'), title, re.IGNORECASE) for x in geneset):
                        g += 1
                    else:
                        continue # if the gene is not in the title, go to next line
                
                # check for unknown species to remove
                if args.unknown is True: 
                    if any(title.upper().startswith(x.upper()) for x in unknown):
                        u += 1 # unknown in title
                        continue # if the unknown is in the title, go to next line
                
                # check the length     
                if minlen is not None and maxlen is not None:    
                    if minlen < int(l) < maxlen:
                        with open(outfile, 'a') as outputfile:
                            print(line, file = outputfile)
                        with open(seqidfile, 'a') as sid:
                            print(seqid, file = sid)
                        s += 1
                elif minlen is not None and maxlen is None:    
                    if minlen < int(l):
                        with open(outfile, 'a') as outputfile:
                            print(line, file = outputfile)
                        with open(seqidfile, 'a') as sid:
                            print(seqid, file = sid)
                        s += 1
                elif minlen is None and maxlen is not None:    
                    if int(l) < maxlen:
                        with open(outfile, 'a') as outputfile:
                            print(line, file = outputfile) 
                        with open(seqidfile, 'a') as sid:
                            print(seqid, file = sid)
                        s += 1
        print(str(t) + ' sequences found')
        if gene is not None:
            print(str(g) + ' sequences with the correct gene found')
        if args.unknown is True: 
            print(str(u) + ' unknown sequences removed')
        print(str(s) + ' sequences passed size selection')
        with open(os.path.join(outfolder, outpfile + '_results.txt'), 'a') as rf:
            print(str(t) + ' sequences found', file=rf)
            print(str(s) + ' sequences passed size selection', file=rf)
        os.remove(os.path.join(outfolder, outpfile + '_selected_taxid.txt'))
        os.rename(os.path.join(outfolder, outpfile + '_sizeselected.txt'), 
                  os.path.join(outfolder, outpfile + '_selected_taxid.txt'))
    else: # no size selection wanted
        with open(os.path.join(outfolder, outpfile + '_selected_taxid.txt'), 'r') as ss:
            for line in ss:
                line = line.strip()
                _, _, _, seqid, title = line.split(' ', 4)
                t += 1 #counting number of sequences scanned
                
                # check if it is the correct gene
                if gene is not None:
                    if any(re.search(x.replace(' ', '.*'), title, re.IGNORECASE) for x in geneset):
                        g += 1
                    else:
                        continue # if the gene is not in the title, go to next line
                
                # check for unknown species to remove
                if args.unknown is True: 
                    if any(title.upper().startswith(x.upper()) for x in unknown):
                        u += 1 # unknown in title
                        continue # if the unknown is in the title, go to next line
                        
                with open(seqidfile, 'a') as sid:
                    print(seqid, file = sid)
        print(str(t) + ' sequences found')
        if gene is not None:
            print(str(g) + ' sequences with the correct gene found')
        if args.unknown is True: 
            print(str(u) + ' unknown sequences removed')
#==============================================================================
# 4. READ REMOVED SEQUENCES FROM A PREVIOUS RUN
# No need to include those again for comparison
#==============================================================================        
def read_removed(p):
    """to read in files from a different pc, rename the files from the other pc
    to something like '0nt_removeset.pick' instead of nt_removeset.pick
    """
    try: # check dbase type
        with open(os.path.join(outfolder, outpfile + '_dbtype.tmp'), 'r') as dbt:
            dbtype = dbt.readline().strip()
    except FileNotFoundError:
        dbtype = check_db_type(outpfile)

    removeset = set()
    try:    # load previous removed sequences, add new ones and store in file
        for dirpath, dirnames, filenames in os.walk('reduced'):
            filenames = [i for i in filenames if i.endswith(dbtype + '_removeset.pick')]
            for name in filenames:
                with open(os.path.join(dirpath, name), 'rb') as f:
                    removeset.update(set(pickle.load(f)))
                os.remove(os.path.join(dirpath, name))
        if p == 'y':
            print('Reading removed accessions from a previous run')  
    except FileNotFoundError:
        pass
    try:
        rf = open(os.path.join(outfolder, outpfile + '_' + dbtype + '_removeset.txt'), 'r')
        if p == 'y':
            print('Reading removed accessions from this run')
        line = rf.readline().strip()
        while line:
            removeset.add(line)        
            line = rf.readline().strip()
        rf.close()
        os.remove(os.path.join(outfolder, outpfile + '_' + dbtype + '_removeset.txt'))
    except FileNotFoundError:
        pass
    if p == 'y':
        print('Saving removed accessions in an updated file')
    outpfilepick = os.path.join(outfolder, outpfile + '_' + dbtype + '_removeset.pick')
    globalpick = os.path.join(outfolder, dbtype + '_removeset.pick')
    with open(outpfilepick, 'wb') as f:
        pickle.dump(removeset, f) # save temp file
    os.replace(outpfilepick, globalpick) # atomic replace of file to global pick file
    return removeset
#==============================================================================
# 5. Read nt or nr fasta file.
# Save reads in species files.  
#==============================================================================
def split_big_files():
    a = ['16S', '16S rRNA', '16S ribosomal RNA',  'small subunit ribosomal', '16S rDNA']
    print('trying to split files bigger than 2 GB to reduce memory consumption') 
    for dirpath, dirnames, filenames in os.walk('reduced'):
        filenames = [i for i in filenames if i.endswith('.spec')]
        for name in filenames:
            if os.path.getsize(os.path.join(dirpath, name))/(1024*1024*1024) > 2: 
                print('processing ' + os.path.join(dirpath, name))
                comparelist = []
                inf = os.path.join(dirpath, name)
                with open(inf, 'r') as infile:
                    for record in SeqIO.parse(infile, "fasta"):
                        comparelist.append([record.id, len(record.seq)])
                comparelist.sort(key=lambda x: x[1]) #sort list based on length seq
                sublist = []
                b = 0
                for i, x in enumerate(comparelist[:-1]):
                    if x[1] > 3000: # reads bigger than 3000 bp, allow 1.5x longer to compare with
                        if x[1]*1.5 < comparelist[i+1][1]:
                            sublist.append([x[0] for x in comparelist[b:i+1]])
                            b = i+1
                    else: # reads smaller than 3000 bp, allow 3x longer to compare with
                        if x[1]*3 < comparelist[i+1][1]:
                            sublist.append([x[0] for x in comparelist[b:i+1]])
                            b = i+1
                sublist.append([x[0] for x in comparelist[b:]])
                
                sublist = [set(x) for x in sublist] 
                
                with open(os.path.join(dirpath, name), 'r') as infile:
                    for record in SeqIO.parse(infile, "fasta"):
                        # search for 16S gene in the name (or other, to be extended)
                        if any(re.search(y.replace(' ', '.*'), record.description, re.IGNORECASE) for y in a):
                            outf = inf.replace('tax_id', 'tax_id_16S_' )
                            with open(outf, 'a') as outputfile:
                                SeqIO.write(record, outputfile, "fasta")
                        else:
                            for i, x in enumerate(sublist):
                                if record.id in x:
                                    outf = inf.replace('tax_id', 'tax_id_' + str(i) + '_' )
                                    with open(outf, 'a') as outputfile:
                                        SeqIO.write(record, outputfile, "fasta")
                os.remove(os.path.join(dirpath, name))
#------------------------------------------------------------------------------
# def download_pickles():
#     try:   
   #         link = urllib.request.urlopen('https://github.com/avierstr/amplicon_sorter'\
   #                                       '/blob/master/amplicon_sorter.py').read()
   #         # find the version-date part of the last version on the webpage
   #         datepart = re.compile(r'(version.*?)(\d{4}-\d{2}-\d{2})(.*version of the script)')
   #         x = datepart.search(str(link))
   #         # the 2nd group of the search is the date
   #         latest_version = x.group(2)
   #         # compare the date of this version with the version on the webpage
   #         if version < latest_version:
   #             version_name = 'amplicon_sorter_' + latest_version + '.py' 
   #             # download latest version
   #             urllib.request.urlopen('https://raw.githubusercontent.com/avierstr/'
   #                                    'amplicon_sorter/master/amplicon_sorter.py')
   #             urllib.request.urlretrieve('https://raw.githubusercontent.com/avierstr'
   #                                        '/amplicon_sorter/master/amplicon_sorter.py',
   #                                        version_name)
   #             print('\n =====================================================\n'
   #                   '| NEW VERSION OF AMPLICON_SORTER AVAILABLE            |\n'
   #                   '| https://github.com/avierstr/amplicon_sorter         |\n'
   #                   '| Downloaded latest version as:                       |\n' 
   #                   '|      ' + version_name + '                  |\n'
   #                   '| Press ctrl-c to exit                                |\n'
   #                   ' =====================================================\n')
   #             t = 10
   #             while t > 0:
   #                 print('Will continue in ' + str(t) + ' seconds...', end='\r')
   #                 time.sleep(1)
   #                 t -= 1
   #             # to clear previous line completely   
   #             print('                                                ', end='\r') 
   #     except:
   #         pass

#------------------------------------------------------------------------------
def sort_species(): # read the new file with exception of the removed accessions
    savefolder = outpfile.replace('.fasta', '')
    try:# check if run was interupted and has to be continued
        inprogressfile = open(os.path.join(outfolder, outpfile + '_inprogress.tmp'), "r") 
        inprogressfile.close()
    except FileNotFoundError: # file not found, so fresh start
        # unknown = open(os.path.join(outfolder, 'unknown_acc_' + filename + '.txt'), 'w') # create empty file
        # unknown.close()
        reducefile = open(os.path.join(outfolder, outpfile + '_reduced.fasta'), 'w') # create empty file
        reducefile.close()        
        
        print('Making temporary folders ready')
        for n in range(0,100): # make temporary folders to store files 
            if n < 10:
                n = '0' + str(n)  #  need folder tmp01 instead of tmp1
            try:
                os.mkdir(os.path.join(outfolder, savefolder, 'tmp'+ str(n))) 
            except OSError:
                #delete subfolder /tmp* and content 
                shutil.rmtree(os.path.join(outfolder, savefolder, 'tmp'+ str(n)))
                time.sleep(0.1)
                os.mkdir(os.path.join(outfolder, savefolder, 'tmp'+ str(n)))
        
        print('----> Reading accessions and taxids')
        taxiddict = {}
        with open(os.path.join(outfolder, outpfile + '_selected_taxid.txt'), 'r') as tid:
            for line in tid:
                tab = line.split() 
                acc, tax_id = tab[0], tab[1] # Assign acc and tax_id  
                taxiddict[acc] = tax_id
                
        removeset = read_removed('y') # check if there are removed accessions from a previous run            
        
        print('Reading fasta file and saving sequences per species without duplicates')
    
        try:
            inputfile = open(os.path.join(outfolder, outpfile + '_selected.fasta'), "r") 
        except FileNotFoundError:
            print('File not found')
        
        t = 0 # number of species in nt or nr file
        for record in SeqIO.parse(inputfile, "fasta"):
            t += 1 #counting number of sequences scanned
            if record.id not in removeset:
                tax_id = taxiddict.get(record.id)
                if tax_id is None: # Genbank uses multiple names if seq are identical
                    d = record.description.split('>')
                    for name in d:
                        # print(name.split()[0])
                        tax_id = taxiddict.get(name.split()[0])
                        if tax_id is not None:
                            record.id = name.split()[0]
                            record.description = name
                            # print(tax_id)
                            break
                foldernr = str(tax_id)[-2:]  #take the last 2 numbers from the taxid
                if len(foldernr) < 2 : # sometimes tax_id is only 1 digit 
                    foldernr = '0' + foldernr  # change to 01
                complete_name = os.path.join(outfolder, savefolder, 'tmp' + foldernr, 'tax_id' + str(tax_id) + '.spec')
                try:
                    with open(complete_name, 'a') as outputfile:
                        SeqIO.write(record, outputfile, "fasta")
                except FileNotFoundError:
                    with open(os.path.join(outfolder, filename + '_unknown_acc.txt'), 'a') as un:
                        print(record.id, file=un)
                    # print(record.description)
                    # print(foldernr)
            if t % 1000000 == 0:  # show progress
                print('\t' + str(t) + ' accessions processed')
                              
        print('Total number of accessions in ' + outpfile + '_selected.fasta: ' + str(t))
        inputfile.close()  
        
        split_big_files() # if there are files bigger than 2 GB, try to split thme in smaller ones
        
        # save a readme file with information about updating a previous reduced file    
        with open(os.path.join(outfolder, 'keep_for_update_README.txt'), 'w') as f:
            print('You can speedup the reduction of the nt or nr database if you keep the', file=f)
            print('files from the last reduction run.', file=f)
            print('- xxx_removeset.pick contains the removed accessions from the previous runs', file=f)
            print('- xxx_uniqueset.pick contains the accessions that have been compared in previous runs', file=f)
            print('Only the newly added accessions will be compared.', file=f)  
         # write a file to indicate work is in progress in case of interruption and restart
        with open(os.path.join(outfolder, outpfile + '_inprogress.tmp'), "w") as inprogressfile:
            print('file to check if analysis is in progress in case of interruption and restart', 
                  file = inprogressfile)
#==============================================================================
# 6. Create a list with official names, TaxID and rank for everything above species level 
#==============================================================================
def check_for_files():
    try: # check if files are downloaded
        with open(os.path.join(dlfolder,'names.dmp'), 'r') as f: 
            pass
    except FileNotFoundError:
        print('Downloading necessary files...')
        ftp = FTP("ftp.ncbi.nlm.nih.gov", timeout=3600) # increase buffersize to 32 MB for ncbi  33554432
        ftp.login("anonymous", "")
        ftp.cwd("/pub/taxonomy/")
        file = 'taxdump.tar.gz'
        filemd5 = file + '.md5'
        print('Downloading ' + file)
        ftp.retrbinary('RETR ' + file , open(os.path.join(dlfolder, file), 'wb').write, blocksize=33554432)
        ftp.retrbinary('RETR ' + filemd5 , open(os.path.join(dlfolder, filemd5), 'wb').write, blocksize=33554432)
        ftp.quit()   
        print('verifying ' + file)
        with open(os.path.join(dlfolder, file), 'rb') as f:
            file_hash = hashlib.md5()
            chunk = f.read(128*1000)
            file_hash.update(chunk)
            while chunk:
                chunk = f.read(128*1000)
                file_hash.update(chunk)
        with open(os.path.join(dlfolder, filemd5), 'r') as f:
            if (f.read().split()[0]) == file_hash.hexdigest():
                print('--> ' + file + ' verified and OK !')
                print('--> Decompressing ' + file )
                my_tar = tarfile.open(os.path.join(outfolder, file), mode='r:gz')  # 'r|gz' data stream of blocks
                my_tar.extract('names.dmp', dlfolder)
                my_tar.extract('nodes.dmp', dlfolder)
                my_tar.close()
                os.remove(os.path.join(dlfolder, file))
                os.remove(os.path.join(dlfolder, file + '.md5'))

#==============================================================================
def read_unique(p):
    with open(os.path.join(outfolder, outpfile + '_dbtype.tmp'), 'r') as dbt:
        dbtype = dbt.readline().strip()
    uniqueset = set()
    try:    # load previous removed sequences, add new ones and store in file
        for dirpath, dirnames, filenames in os.walk('reduced'):
            filenames = [i for i in filenames if i.endswith(dbtype + '_uniqueset.pick')]
            for name in filenames:
                with open(os.path.join(dirpath, name), 'rb') as f:
                    uniqueset.update(set(pickle.load(f)))
                os.remove(os.path.join(dirpath, name))
        if p == 'y':    
            print('Reading unique accessions from a previous run')  
    except FileNotFoundError:
        pass
    try:
        rf = open(os.path.join(outfolder, outpfile + '_' + dbtype + '_uniqueset.txt'), 'r')
        if p == 'y':
            print('Reading unique accessions from this run')
        line = rf.readline().strip()
        while line:
            uniqueset.add(line)        
            line = rf.readline().strip()
        rf.close()
        os.remove(os.path.join(outfolder, outpfile + '_' + dbtype + '_uniqueset.txt'))
    except FileNotFoundError:
        pass
    if p == 'y':
        print('Saving unique accessions in an updated file')
    outpfilepick = os.path.join(outfolder, outpfile + '_' + dbtype + '_uniqueset.pick')
    globalpick = os.path.join(outfolder, dbtype + '_uniqueset.pick')
    with open(outpfilepick, 'wb') as f:
        pickle.dump(uniqueset, f) # save temp file
    os.replace(outpfilepick, globalpick) # atomic replace of file
    return uniqueset, dbtype
#==============================================================================
def cleanup():
    cl = True
    while cl is True:
        time.sleep(30)
        tmpdirs = []
        for dirpath, dirnames, filenames in os.walk(outfolder):
            for name in dirnames:
                if name.startswith('tmp'):
                    tmpdirs.append((dirpath, name))
        for p, n in tmpdirs:
            pn = os.path.join(p, n)
            if not os.listdir(pn):
                os.rmdir(pn)
        if len(tmpdirs) == 0:
            cl = False
#==============================================================================    
def process_file(todoqueue, dbtype):  # read each species file and compare sequences
    for infile in iter(todoqueue.get, 'STOP'):    #do stuff until infile.get returns "STOP"
        try:      
            starttime = time.time()
            MYLOCK = Lock()  
            t = 0 #number of sequences
            s = 0 #number of similarities found  
            z = 0 #number in list to show progress
            comparelist = [] # list for comparing the sequences
            uniquelist = [] # list with unique sequences after comparison
            removelist = [] # list with sequences removed after comparison
            uniquesetname = os.path.join(outfolder, outpfile + '_' + dbtype + '_uniqueset.txt')
            removesetname = os.path.join(outfolder, outpfile + '_' + dbtype + '_removeset.txt')
            inputfile = open(infile, "r") 
            for record in SeqIO.parse(inputfile, "fasta"):
                t += 1 #counting number of sequences scanned   
                if record.id in uniqueset:
#                    print('in unique')
                    comparelist.append([record.id, str(record.seq), 'u']) #unique
                elif record.id in removeset:  # if acc was removed in interrupted run
#                    print('in removed')
                    pass
                    # comparelist.append([record.id, ' ', ''])  # no need to compare again, need it for position to restart
                else:
#                    print('new')
                    comparelist.append([record.id, str(record.seq), 'n']) #new
                    
            comparelist.sort(key=lambda x: len(x[1])) #sort list based on length seq
            inputfile.close()
            
            folder, name = os.path.split(infile)
            lastf = folder.split('/')[-1]
            f = '/'.join(('', lastf, name))
            
            print('Processing ' + f + ' containing ' + str(t) + ' sequences.') 
            # compare 1 with 2,3,4,5,...; compare 2 with 3,4,5,... compare 3 with 4,5,...

            # temp file to know untill where comparisons have been done if run was interrupted
            comparisonfile = infile[:-5] + '.comp'  
            try:   # if interrupted, start where left off
                cf = open(comparisonfile, 'r')
                comparison = cf.readline()
                cf.close()
                C1 = comparison.strip()
                try:
                    position = ([n for n, [acc, seq, val] in enumerate(comparelist) if acc == C1])[0]
                    z = B1 = position
                    print('restarting at position ' + str(position))
                except IndexError:
                    B1 = ''        #initialize values for comparisonfile  
                    position = 0
            except FileNotFoundError:
                B1 = ''        #initialize values for comparisonfile  
                position = 0  #start from first item in list
                        
            for position in range(position, len(comparelist)-1):
                z += 1
                if time.time()-starttime > 1800:  # show progress
                    starttime = time.time()
                    print('---> ' + infile[7:] + ' (' + str(t) + ' seq) ' +
                          str(round(z/t*100)) + '% done... (' + str(s) + ' similar)')
                    # store temporary results in case it takes days to process file 
                    # only removed sequences are sure
                    removelist = [i[0] for i in comparelist if i[1] == ' ']
                    MYLOCK.acquire()
                    with open(removesetname, 'a') as rf:
                        for x in removelist:
                            print(x, file=rf)
                    # store the last comparison in a file
                    with open(comparisonfile, 'w') as cf:
                        print(B1, file=cf)
                    # all sequences before "position" that are not removed are unique
                    uniquelist = [i[0] for i in comparelist[:position] if i[1] != ' ']
                    with open(uniquesetname, 'a') as uf:
                        for x in uniquelist:
                            print(x, file=uf)
                    MYLOCK.release()
                    
                for position2 in range(position+1,len(comparelist)):
                    A1 = comparelist[position][1]
                    A2 = comparelist[position2][1]
                    B1 = comparelist[position][0]
#                    print(comparelist[position][2] + ' - ' + comparelist[position2][2])
#                    print(str(len(A1)) + ' - ' + str(len(A2)))
                    if comparelist[position][2] == 'u' and comparelist[position2][2] == 'u':  
#                        print('both were compared')
                        break    #compared in previous run, no need to do it again
                    elif A1 == ' ': 
#                        print('A1 was removed')
                        break  # if one of the sequences is removed, no need to compare
                    elif len(A1)*3 < len(A2):
                        break # don't compare is length difference is too big
                    else:  
                        iden = distance(A1, A2)
                        if iden > 0.97: 
                            # print(str(len(A2)/len(A1)) + ' - ' + str(len(A1)) + ' - ' + str(len(A2)))
                            # print(comparelist[position][0] + ' - ' + comparelist[position2][0])
                            s += 1                           
                            comparelist[position][1] = ' ' 

            print('File ' + f +': ' + str(s) + '/' + str(t) + ' similar sequences found')
            
            # only keep unique and removed acc numbers
            uniquelist = [i[0] for i in comparelist if i[1] != ' ']
            removelist = [i[0] for i in comparelist if i[1] == ' ']
            
            MYLOCK.acquire()
            with open(os.path.join(outfolder, 'reduce.log'), 'a') as lf:
                print('File ' + f +': ' + str(s) + '/' + str(t) + ' similar sequences found', file=lf)
                
            with open(uniquesetname, 'a') as uf:
                 for x in uniquelist:
                     print(x, file=uf)
            
            with open(removesetname, 'a') as rf:
                 for x in removelist:
                     print(x, file=rf)
                     
            rd = 0 #number of reduces sequences
            inputfile = open(infile, "r")
            reducedfilename = os.path.join(outfolder, outpfile + '_reduced.fasta')
            for record in SeqIO.parse(inputfile, "fasta"):
                if record.id in uniquelist:
                    rd += 1          
                    with open(reducedfilename, 'a') as outputfile:
                        SeqIO.write(record, outputfile, "fasta")
            inputfile.close()

#            print('process_id ' + str(b.pid) +' ' + str(rd) + ' unique sequences added to ' + filename + '_red.fasta')
            os.remove(infile)
#            print('process_id ' + str(b.pid) + ' --> ' + infile + ' deleted')

            try:  # remove temporary file if exists
                os.remove(comparisonfile)
            except FileNotFoundError:
                pass
            
            MYLOCK.release()

        except KeyboardInterrupt:
            print("Shutting processes down")
            sys.exit()
#==============================================================================
def compare(dbtype): # process all the files with multiprocessing
    nprocesses = args.nproc
    print('Warming up, meanwhile can you get me a coffee ?')
    print('Starting ' + str(nprocesses) + ' processes')
    todoqueue = Queue()#maxsize = nprocesses*4) # keep the queue max 4 times the number of processes
    def queuer(): # feed the queue
        for dirpath, dirnames, filenames in os.walk(outfolder):
            dirnames.sort()
            filenames = [i for i in filenames if i.endswith('.spec')]
            for name in filenames:
                todoqueue.put(os.path.join(dirpath, name))#, block=True)
        for i in range(nprocesses): # put 'STOP' at the end of the queue for every process
                todoqueue.put("STOP")   
    def consumer():
        try:
            process = [Process(target=process_file, args=(todoqueue, dbtype,)) for x in range(nprocesses)]
#            for p in process:
                #ask the processes to stop when all files are handled
                #"STOP" is at the very end of queue
#                todoqueue.put("STOP")
            for p in process:
                p.start()
            for p in process:
                p.join()    
        except KeyboardInterrupt:                
            print("Shutting processes down")
            sys.exit()
           # Optionally try to gracefully shut down the worker processes here.       
#            p.terminate()
#            p.join()          
    cl = Thread(target = cleanup)
    cl.start()
    c = Thread(target = consumer)
    c.start()
    q = Thread(target = queuer)
    q.start()
    q.join()
    c.join() # wait until c has finished its work
    cl.join()
    
    inputfile = open(os.path.join(outfolder, outpfile + '_reduced.fasta'), "r") 
    t = 0
    for record in SeqIO.parse(inputfile, "fasta"):
        t += 1 #counting number of sequences
    inputfile.close()
    print(str(t) + ' sequences in the reduced ' + outpfile + ' database')    
    with open(os.path.join(outfolder, outpfile + '_results.txt'), 'a') as rf:
        print(str(t) + ' sequences in the reduced ' + outpfile + ' database', file=rf) 
    os.remove(os.path.join(outfolder, outpfile + '_inprogress.tmp'))
#==============================================================================
# 7. Create BLAST database from reduced file
#==============================================================================
def make_blast_db(outpfile):
         
    try: # check if reduced file exits
        f = open(os.path.join(outfolder, outpfile + '_reduced.fasta'), 'r')
        f.close()
        inf = outpfile + '_reduced.fasta'
    except FileNotFoundError:
        inf = outpfile + '_selected.fasta'
    print('Creating BLAST database from ' + inf)
    # create taxid_map
    taxid_mapf = os.path.join(outfolder, outpfile + '_taxid_map.txt')
    f = open(taxid_mapf, 'w')
    with open(os.path.join(outfolder, outpfile + '_selected_taxid.txt'), 'r') as ss:
        for line in ss:
            line = line.strip()
            acc, tax, l, seqid, title = line.split(' ', 4)
            # with open(taxid_mapf, 'a') as tm:
                # print(seqid + ' ' + tax, file = tm)
            f.write(seqid + ' ' + tax + '\n')
    f.close()
            
    savefolder = outpfile.replace('.fasta', '')
    try: # check dbase type
        with open(os.path.join(outfolder, outpfile + '_dbtype.tmp'), 'r') as dbt:
            dbtype = dbt.readline().strip()
    except FileNotFoundError:
        dbtype = check_db_type(outpfile)
        
    command = ['makeblastdb']
    command.append('-in')
    command.append(os.path.join(outfolder, inf))
    command.append('-title')
    command.append(inf)
    command.append('-parse_seqids')
    command.append('-taxid_map')
    command.append(taxid_mapf)
    command.append('-out') 
    command.append(os.path.join(outfolder, savefolder, outpfile))
    command.append('-dbtype')
    command.append(dbtype)
    command.append('-max_file_sz')
    command.append('3GB')
    command.append('-hash_index')
    # print(command)
    result = subprocess.run(command, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0 :
        print('----> Error ' + str(result.stderr) + '\n')
    else:
        print('----> Blast database created !')  
        
    os.remove(os.path.join(outfolder, outpfile + '_dbtype.tmp'))
#==============================================================================
if __name__ == '__main__':
    check_version(version)
    args = arguments()

    if args.subparser_name == 'tax':
        folder = args.outfolder
        # create an outputfolder
        if folder is not None:
            outfolder = folder
        else:
            outfolder = os.getcwd()
        if not os.path.exists(outfolder):
            os.mkdir(outfolder) # create the folder
            
        if args.taxidlist is not None: # create a taxidlist for reference
            create_taxidlist(args.taxidlist)
        if args.taxtree is not None: # create a taxidtree for reference
            create_taxtree(args.taxtree)
        if args.lineage is not None: # create a lineage for reference
            create_lineage(args.lineage)
    else:
        folder = args.outfolder
        # create an outputfolder for download
        if folder is not None:
            dlfolder = folder
            if not os.path.exists(dlfolder):
                os.mkdir(dlfolder) # create the folder
        else:
            dlfolder = os.getcwd()
            
        # create an outputfolder for sub-databases
        if folder is not None:
            outfolder = os.path.join(folder, 'reduced')
        else:
            outfolder = os.path.join(os.getcwd(), 'reduced')
        if not os.path.exists(outfolder):
            os.mkdir(outfolder) # create the folder
            
        filename = args.database
        
        filenames = os.listdir(outfolder)
        filenames = [n for n in filenames if n.endswith('_inprogress.tmp')]
        if len(filenames) > 0 or args.append is True or args.select is not None \
            or args.reduce is True or args.makeblastdb is True:
            outpfile = args.outfile
            if outpfile is None:
                print('Please provide the output filename with "-o" or "-outfile"')
                sys.exit()
            savefolder = os.path.join(outfolder, outpfile.replace('.fasta', ''))
            if not os.path.exists(savefolder):
                os.mkdir(savefolder) # create the folder
            save_arguments() # write all settings in the results.txt file
            
        if len(filenames) == 0: # new run
            if args.download is True:  # download necessary files
                get_files(filename)
            
            if args.append is True:
                original_outpfile = outpfile
                outpfile = outpfile + '_append'
                
            if args.select is not None: # select specific species
                taxids = args.select
                select_sequences(taxids, outpfile)
                if args.append is True:
                    with open(os.path.join(outfolder, original_outpfile + '_results.txt'), 'a') as rf:
                        with open(os.path.join(outfolder, outpfile + '_results.txt'), 'r') as infile:
                            content = infile.read()
                            rf.write(content)
                    with open(os.path.join(outfolder, original_outpfile + '_selected.fasta'), 'a') as rf:    
                        with open(os.path.join(outfolder, outpfile + '_selected.fasta'), 'r') as infile:
                            for line in infile:
                                rf.write(line)
                    with open(os.path.join(outfolder, original_outpfile + '_selected_taxid.txt'), 'a') as rf:    
                        with open(os.path.join(outfolder, outpfile + '_selected_taxid.txt'), 'r') as infile:
                            for line in infile:
                                rf.write(line)
                    os.remove(os.path.join(outfolder, outpfile + '_dbtype.tmp'))
                    os.remove(os.path.join(outfolder, outpfile + '_results.txt'))
                    os.remove(os.path.join(outfolder, outpfile + '_selected.fasta'))
                    os.remove(os.path.join(outfolder, outpfile + '_selected_taxid.txt'))
                    outpfile = original_outpfile
                    
            if args.reduce is True: # remove highly similar sequences from each species
                # download_pickles()
                sort_species()  # read the new file with exception of the removed accessions
                # start comparing sequences.  Read remove and unique again in case run was interrupted
                removeset = read_removed('y') # check if there are removed accessions from a interrupted run
                uniqueset, dbtype = read_unique('y') # check if there are unique accessions from a previous run
                compare(dbtype)
                print('Cleaning up...')
                _ = read_removed('n') # check if there are removed accessions from a interrupted run
                _ = read_unique('n') # check if there are unique accessions from a previous run
                print('Done.') 
        else: # run was interrupted
            print('Interrupted "reduce" run detected, continue reducing')
            for progressname in filenames:
                if outpfile == progressname.replace('_inprogress.tmp', ''):
                    removeset = read_removed('y') # check if there are removed accessions from a interrupted run
                    uniqueset, dbtype = read_unique('y') # check if there are unique accessions from a previous run
                    compare(dbtype)
                    print('Cleaning up...')
                    _ = read_removed('n') # check if there are removed accessions from a interrupted run
                    _ = read_unique('n') # check if there are unique accessions from a previous run
                    print('Done.') 
        
        if args.makeblastdb is True:
            make_blast_db(outpfile)
            