#!/usr/bin/env python
from Bio import SeqIO
import sys
import csv
import argparse

parser = argparse.ArgumentParser(description='''The script takes a matrix of asvs with the asvs as the rows and the samples in columns''')
parser.add_argument('-tax_ref', help = 'the path to the  reference taxonomy table that matches the ref fasta file searched agaist with a taxonomy string e.g. silva.tax')
parser.add_argument('-hits', help = 'the output from vsearch when run like this - vsearch --usearch_global dada2-fasta.fa --db /groups/g454/blastdbs/gast_distributions/silva119.fa --blast6out NODE-HITS.txt --id 0.6  *** keep in mind that the --db is dependent on the gene of interest')
parser.add_argument('-dada2', help = 'the count or percent abundance matrix produced by Dada2 and edited asv names as rows and samples in columns the A1 position must be "samples"')
parser.add_argument('-fa', help = 'the fasta of dada2 representative sequences')
args = parser.parse_args()

def open_node_hits_to_dict(sample_name):
    sample_dict = {}
    with open(sample_name, 'r') as f:
        for line in f:
            x = line.strip().split("\t")
            sample_dict[x[0]] = x[1]
    return sample_dict

def open_silva_database_to_dict(db_name):
    db_name_dict = {}
    with open(db_name, 'r') as f:
        for line in f:
            x = line.strip().split('\t')
            db_name_dict[x[0]] = x[1:len(x)]
    return db_name_dict

def open_dada2_table_to_dict(med_table):
    med = {}
    with open(med_table, 'r') as f:
        for line in f:
            x = line.strip().split('\t')
 #           print(x)
            med[x[0]] = x[1:len(x)]
    return med

def open_NODE_REPRESENTATIVES(node_fasta):
    all_nodes = {}
    with open(node_fasta, 'r') as f:
        x = SeqIO.parse(f, "fasta")
        for rec in x:
            all_nodes[rec.id] = [rec.id]
    return all_nodes

# Create each of the dictionaries using the defined functions above
node_hits_dict = open_node_hits_to_dict(args.hits)
tax_lookup_dict = open_silva_database_to_dict(args.tax_ref)
dada2_matrix_dict = open_dada2_table_to_dict(args.dada2)
all_nodes_dict = open_NODE_REPRESENTATIVES(args.fa)

#for key in dada2_matrix_dict.keys():
#    print(key,dada2_matrix_dict[key])

#Open and output file
output_tax = open('PHYLOSEQ-TAX.txt', 'w')
output_matrix = open('PHYLOSEQ-MATRIX.txt', 'w')

#Add header to output file
output_tax.write("node"+';'+"db_hit_id"+';'+"Kingdom"+';'+"Phylum"+';'+"Class"+';'+"Family"+';'+"Order"+';'+"Genus"+';'+"Species"+";"+"strain"+";"+"strain1"+";"+"ASV"+'\n')

# Add the node and taxonomy string to the output file
# Add the nodes that rec'd no hit in the  silva db to the taxonomic output.
for key in all_nodes_dict.keys():
    if key in node_hits_dict.keys():
         output_tax.write(str(key)+';'+str(node_hits_dict[key])+';'+str(tax_lookup_dict[node_hits_dict[key]][0])+';'+str(key)+'\n')
    else:
        output_tax.write(str(key)+';'+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+"na;"+'\n')

# Write the header information to the new matrix file
output_matrix.write("node"+'\t'+'\t'.join(dada2_matrix_dict['sample'])+'\n')

# Write the transposed phyloseq matrix to the output file
for key in all_nodes_dict.keys():
    if key in dada2_matrix_dict.keys():
#        print(key,dada2_matrix_dict[key], "wow")
        output_matrix.write(str(key)+'\t'+'\t'.join(dada2_matrix_dict[key])+'\n')
output_tax.close()
output_matrix.close()
