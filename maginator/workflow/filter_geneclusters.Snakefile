#!/bin/bash

from itertools import product
import os

WD = config['wd']
F_DIR = os.path.join(WD, 'clusters')
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

CLUSTERS = set(glob_wildcards(os.path.join(F_DIR, '{cluster}')).cluster)
CLUSTERS = {x for x in CLUSTERS if x.isdigit()}

SAMPLES = set()
FASTQ1 = set()
FASTQ2 = set()
sample_dict = {}

with open(param_dict['reads']) as f:
    for row in f:
        row = row.rstrip().split(',')
        SAMPLES.add(row[0])
        sample_dict[row[0]] = [row[1], row[2]]
        FASTQ1.add(row[1])
        FASTQ2.add(row[2])

# make dict of lists of authorized wildcard combinations (based on presence of input files) - dict with key as clusters with value as list of samples
my_comb={}
for comb in product(SAMPLES, CLUSTERS):
    path = os.path.join(WD, 'clusters')
    
    # going through all the files with combinations of clusters and samples saving it in a dict
    if(os.path.isfile(path + "/" + comb[1] + "/" + comb[0] + "_" +comb[1] + ".fa")):
        if comb[1] in my_comb: # if the cluster is present in the dict add the extra sample
            my_comb[comb[1]].append(comb[0])
        else: # if the cluster is new add the cample and cluster to the dict
            my_comb[comb[1]] = [comb[0]]


rule all:
    input:
        os.path.join(WD, 'genes', 'redundant_genes.tsv'),
        expand(os.path.join(WD, 'clusters', 'gene_lists', '{cluster}_genes.txt'), cluster=CLUSTERS),
        expand(os.path.join(WD, 'mapped_reads', 'gene_counts_{sample}.bam'), sample=SAMPLES)


# creating a file only containing the redundant genes
# Based on the clustering file from mmseqs
rule redundant_genes:
    input:
        os.path.join(WD, 'genes', 'all_genes95_cluster.tsv')
    output:
        os.path.join(WD, 'genes', 'redundant_genes.tsv')
    resources:
        cores = 1,
        memory = 32,
        runtime = '1:00:00'
    shell: 
        "awk '$1!=$2' {input} | cut -f2 | uniq > {output}"


# Creating a list with the nonredundant genes for each cluster only containing non-redundant genes (also excluding representatives) 
rule nonredundant_clusters:
    input:
        nuc = lambda wildcards: expand(os.path.join(WD, 'clusters', '{{cluster}}', '{sample}_{{cluster}}.fa'), sample=my_comb[wildcards.cluster]),
        redundant = os.path.join(WD, 'genes', 'redundant_genes.tsv')
    output:
        nuc_nonredundant = temp(os.path.join(WD, 'clusters', 'nonredundant', '{cluster}_nonredundant.fna')),
        gene_lists = os.path.join(WD, 'clusters', 'gene_lists', '{cluster}_genes.txt')
    resources:
        cores = 1,
        mem_gb = 50,
        runtime = 86400
    shell:
        "cat {input.nuc} | filterbyname.sh in=stdin.fa out={output.nuc_nonredundant} names={input.redundant} include=f; "
        "grep '>' {output.nuc_nonredundant} | cut -d' ' -f1 | sed --expression='s/>//g' > {output.gene_lists} || true "

# Creating a nonredundant catalogue of all the genes
rule nonredundant_catalogue:
    input:
        genecat = os.path.join(WD, 'genes', 'all_genes.fna'),
        redundant = os.path.join(WD, 'genes', 'redundant_genes.tsv')
    output:
        genecat_nonredundant = os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 50,
        runtime = '12:00:00'
    shell:
        "filterbyname.sh in={input.genecat} out={output.genecat_nonredundant} names={input.redundant} include=f"

# Indexing the genes for mapping
rule bwa_index:
    input: os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta')
    output: done =touch(os.path.join(WD, 'genes', 'all_genes_nonredundant')),
            gene_lengths = os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta.fai')
    conda:
        "envs/filter_geneclusters.yaml" 
    resources:
        cores = 40,
        mem_gb = 188, 
        runtime = '1:00:00:00'
    shell:
        "bwa-mem2 index {input}; samtools faidx {input}; touch {output.gene_lengths}"

# Readmapping
rule bwa_readmap:
    input:
        index = os.path.join(WD, 'genes', 'all_genes_nonredundant'),
        gene_cat = os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta'),
        fastq1 = lambda wildcards: sample_dict[wildcards.sample][0],
        fastq2 = lambda wildcards: sample_dict[wildcards.sample][1]
    params:
        sample = SAMPLES
    output:
        bam = os.path.join(WD, 'mapped_reads', 'gene_counts_{sample}.bam')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 40,
        mem_gb = 188,
        runtime = '1:00:00:00'
    shell:
        "bwa-mem2 mem -t {resources.cores} {input.gene_cat} {input[2]} {input[3]} | samtools view -T {input.gene_cat} -F 3584  -b --threads {resources.cores} | samtools sort --threads {resources.cores} > {output.bam};"
