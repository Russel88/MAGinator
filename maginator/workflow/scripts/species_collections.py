#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2022

@author: trinezachariasen
"""

# importing the species collections and it's corresponding MGS-ids included in the Metagenomic Species Collection
import csv
import sys
csv.field_size_limit(sys.maxsize) # I have large entries in my csv-file with gene-id's

cluster_MGS_d = {}
with open(snakemake.input[0] , 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
        for cluster in row[2].split(','):
            cluster_MGS_d[cluster] = row[0]

d = {}
with open(snakemake.input[0], 'r') as csv_file:
    for row in csv.reader(csv_file, delimiter='\t'):
            d[row[0]] = row[2] # setting the MGS as key and it's included MGS-ids as value

# importing the VAMB clustefile as dict with the geneid as key and it's corresonding MGS as value
id_dict = {}
with open(snakemake.input[2], 'r') as id_file:
    for row in csv.reader(id_file, delimiter='\t'):
        cluster = row[0].split("_")[-1]
        id_dict[row[1]] = cluster

## Creating a file only containing the representative / nonredundant genes
## Based on the clustering file from mmseqs
gene_file = snakemake.input[1]
outfile = snakemake.output[0]


excluded_counter = set()
cluster_counter = set()
contig_gene = {}

with open(gene_file) as f:
    for line in f.readlines():
        [raw_rep, raw_gene] = line.rstrip().split('\t') # the geneIDs
        rep = "_".join(raw_rep.split("_")[:-1]) # representative contigID
        gene = "_".join(raw_gene.split("_")[:-1]) # clustered contigID

        repMGS = id_dict[rep] # mgsID of representative
        geneMGS = id_dict[gene] # mgsID of the gene
        if (repMGS not in cluster_MGS_d):
            cluster_MGS_d[repMGS]=repMGS

        if (rep==gene): #if the representative gene is same as clustered we keep it
            cluster_counter.add(raw_rep)
            contig_gene[raw_rep] = cluster_MGS_d[repMGS]
        try:
            if ((d[geneMGS].find(repMGS)==0)): # if the MGS of the gene is included in the representative gene Metagenomic Species Collection
                cluster_counter.add(raw_rep)
                contig_gene[raw_rep] = cluster_MGS_d[repMGS]
        except KeyError:
            pass

        if (geneMGS==repMGS):  # if the MGS of the Gene is same as representative gene MGS
             cluster_counter.add(raw_rep)
             contig_gene[raw_rep] = cluster_MGS_d[repMGS]
        else: # if the representative MGSid is not equal to the MGS of the gene it represents
            excluded_counter.add(raw_rep) 
            contig_gene.pop(raw_rep, None)
             
# identifying the representative genes, which do not cover more than one cluster
rep_to_keep = list(cluster_counter.difference(excluded_counter))

with open(outfile, "w") as out:
    out.write("\n".join(rep_to_keep))

# sorting the geneID_clusterID to become geneID_collectionID
with open(snakemake.output[1],'w') as f:
    w = csv.writer(f, delimiter='\t')
    for gene in rep_to_keep:
       w.writerow([gene, contig_gene[gene]])

