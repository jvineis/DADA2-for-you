#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(description='''Use the dada2 output edited with the column header "sample", ASVs as rownames and samples as column headers as input to create a separate fasta file and count table with a new asv identifier in the table and a beautiful matching fasta file''')
parser.add_argument('-n', help = 'file containig DNA sequences, each sequence on a new line')
parser.add_argument('-fa', help = 'filename for the beautiful fasta file')
parser.add_argument('-o', help = 'filename for the matching DADA2 count table that uses the new asv ids instead of the sequence as the identifier')
args = parser.parse_args()

infile = open(args.n, 'r')
outtable = open(args.o, 'w')
outfasta = open(args.fa, 'w')


table_dict = {}
for line in infile:
    x = line.strip().split('\t')
    table_dict[x[0]] = x[1:len(x)]

    
for key in table_dict.keys():
    if str(key) == "sample":
        outtable.write("sample"+'\t'+'\t'.join(table_dict[key])+'\n')
    
count = 1
for key in table_dict.keys():
    if key == "sample":
        next
    else:
        outfasta.write(">asv_"+str(count)+'\n'+str(key)+'\n')
        outtable.write("asv_"+str(count)+'\t'+'\t'.join(table_dict[key])+'\n')
    count += 1
