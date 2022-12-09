#!/usr/bin/env python3

import re
import sys

import pandas as pd
import numpy as np

from scipy.sparse import lil_matrix
from Bio import SeqIO

# Get clustering table
dat = pd.read_csv(snakemake.input['cluster'], header=None, names=('Cluster', 'Gene'), sep='\t')

# Get cluster index
cluster_index = list(set(dat['Cluster']))
n_clust = len(cluster_index)

# Turn into dicts
clust_gene = dict(zip(dat['Gene'], dat['Cluster']))
clust_index = dict(zip(cluster_index, range(n_clust)))

# Init graph as sparse array
G = lil_matrix((n_clust, n_clust))

# Read genes and save in graph
k = 0
prev_prog = 0
prev_contig = None
prev_cluster = None
prev_strand = None
n_gene = dat.shape[0]
with open(snakemake.input['faa']) as fh:
    for gene in SeqIO.parse(fh, 'fasta'):
        # Get contig name
        contig = re.sub('_[0-9]*$', '', str(gene.id))
        
        # Get gene cluster
        cluster = clust_gene[str(gene.id)]
       
        # Get strand
        strand = str(gene.description).split()[6]

        # If we are on the same contig and strand as previous gene, genes are adjacent
        if contig == prev_contig and strand == prev_strand:
            # Get indices
            i = clust_index[cluster]
            j = clust_index[prev_cluster]
            G[i, j] += 1
            G[j, i] += 1

        # Save previous cluster
        prev_cluster = cluster

        # Save previous contig
        prev_contig = contig

        # Save previous strand
        prev_strand = strand

        prog = round(100*k/n_gene)
        if prog > prev_prog:
            print(str(prog)+'%')
        prev_prog = prog
        k += 1
      
# Write graph
with open(snakemake.output['graph'], 'w') as fh:
    print(G, file = fh)

# Write cluster indices
with open(snakemake.output['index'], 'w') as fh:
    for k in zip(range(n_clust), cluster_index):
        fh.write('{} {}\n'.format(k[0], k[1]))
