#!/usr/bin/env python3

import glob
import re
import sys
import os

import multiprocess as mp

from Bio import SeqIO

# Constants
align_dir = snakemake.params['ali_dir']
samp_tab = snakemake.input[0]
out_dir = snakemake.output[0]
threads = snakemake.resources['cores']

# Make out dir
try:
    os.mkdir(out_dir)
except FileExistsError:
    pass

# Get list of clusters
clusters = [x for x in os.listdir(align_dir)]

# Get list of samples
samples = {}
samples['Outgroup'] = ''
with open(samp_tab, 'r') as handle:
    for x in handle.readlines():
        samples[x.strip().split(',')[0]] = ''

# For each cluster:
'''
Load sequences
Concat
Write to a single fasta file
Write partition file for iqtree 
'''

def concat(clust, samples = samples):

    # Wrapper to write partition file
    with open(os.path.join(out_dir, clust+'.part'), 'w') as part_handle:

        # Init
        gene_start = 1
        seq_len = 1

        # Read through alignments and fill in
        for afile in glob.glob(os.path.join(align_dir, clust, clust+'_'+'*.fna')):
            seqs = {}

            # Load in seqs
            with open(afile, 'r') as handle:
                for fa in SeqIO.parse(handle, 'fasta'):
                    seqs[fa.id] = str(fa.seq)
                    seq_len = len(fa.seq)
           
            # Loop through samples and fill in
            for samp in samples:
                if samp in seqs.keys():
                    samples[samp] += seqs[samp]
                else:
                    samples[samp] += '-' * seq_len

            # Write partition file
            gene_name = re.sub('.*'+clust+'_(.*)\.fna$', r'\1', afile)
            gene_end = gene_start + seq_len - 1
            part_handle.write('DNA, '+str(gene_name)+' = '+str(gene_start)+'-'+str(gene_end)+'\n')
            gene_start += seq_len

    # Write fasta
    samples = {k:v for (k,v) in samples.items() if len(v.replace('-', '')) > 0}

    with open(os.path.join(out_dir, clust+'.fna'), 'w') as handle:
        for samp in samples:
            handle.write('>'+samp+'\n'+samples[samp]+'\n')

pool = mp.Pool(threads)
concats = list(pool.imap(concat, set(clusters)))

