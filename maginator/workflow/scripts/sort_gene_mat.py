#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 09:59:31 2022

@author: trinezachariasen
"""

import csv
import sys
csv.field_size_limit(sys.maxsize) # I have large entries in my csv-file with gene-id's


# importing set with the genes we want to keep
gene_set = set()
with open(snakemake.input[1], 'r') as id_file:
    for row in csv.reader(id_file):
        gene_set.add(row[0])
    
out_mat = open(snakemake.output[0],"w")

# Going through the gene count matrix (containing all genes and all clusters, sorting away the genes, that span multiple clusters)
with open(snakemake.input[0], 'r') as mat_file:
    out_mat.write(mat_file.readline()) # writing the header to the outfile
    
    for row in csv.reader(mat_file, delimiter='\t'): # looping through the genes
#        gene = "_".join(row[0].split("_")[:-1])
        if (row[0] in gene_set): #if the gene is found in the gene_set we want to keep, write it to the smaller matrix-file
            out_mat.write("\t".join(row) + "\n")

