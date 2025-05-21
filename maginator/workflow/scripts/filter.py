import os
import sys
import re
import logging

from Bio import SeqIO

'''
Filter bins by size
1st input VAMB clusters.tsv file
2nd input contigs fasta file
3rd input output directory for writing a fasta per cluster
4th input is bin size limit
'''

# Read VAMB file into dict
with open(snakemake.input[0], 'r') as fh:
    vamb_dict = [x.strip().split() for x in fh.readlines()]
vamb_dict = {x[1]:x[0] for x in vamb_dict}

# Read contigs and organize into dict per bin
bin_contigs = dict()
with open(snakemake.input[1], 'r') as fh:
    for fa in SeqIO.parse(fh, 'fasta'):
        bin_id = vamb_dict[fa.id]
        try:
            bin_contigs[bin_id].append(fa)
        except Exception:
            bin_contigs[bin_id] = [fa]

# Filter by length
bin_len = {x: sum([len(y.seq) for y in bin_contigs[x]]) for x in bin_contigs}
bin_filter = [x for x in bin_len if bin_len[x] >= int(snakemake.params.binsize)]

# Write per cluster
cluster_set = set([int(re.sub('.*_', '', x)) for x in bin_filter])

try:
    os.mkdir(snakemake.output[0])
except Exception:
    pass

for cluster in cluster_set:
    bins = [x for x in bin_filter if bool(re.search('_'+str(cluster)+'$', x))]
    clust_path = os.path.join(snakemake.output[0], str(cluster))

    try:
        os.mkdir(clust_path)
        os.mkdir(snakemake.output[1])
    except Exception:
        pass

    for bin_id in bins:
        SeqIO.write(bin_contigs[bin_id], os.path.join(snakemake.output[0], str(cluster), bin_id+'.fa'), "fasta")
        SeqIO.write(bin_contigs[bin_id], os.path.join(snakemake.output[1], bin_id+'.fa'), "fasta")

