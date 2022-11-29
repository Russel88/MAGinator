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
sample_dict = {}

with open(param_dict['reads']) as f:
    for row in f:
        row = row.rstrip().split(',')
        SAMPLES.add(row[0])
        sample_dict[row[0]] = [row[1], row[2]]

rule all:
    input:
        expand(os.path.join(WD, 'mapped_reads', 'bams', 'gene_counts_{sample}.bam'), sample=SAMPLES)

# Creating a nonredundant catalogue of all the genes
rule nonredundant_catalogue:
    input:
        genecat = os.path.join(WD, 'genes', 'all_genes.fna'),
        clusters = os.path.join(WD, 'genes', 'all_genes_cluster.tsv')
    output:
        os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta')
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores = 1,
        memory = 50,
        runtime = '12:00:00'
    shell:
        "perl -ne 'if(/^>(\S+)/){{$c=$i{{$1}}}}$c?print:chomp;$i{{$_}}=1 if @ARGV' <(cut -f1 {input.clusters} | uniq) {input.genecat} | awk '{{print $1}}' > {output}"

# Indexing the genes for mapping
rule bwa_index:
    input: 
        fasta=os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta')
    output: 
        index=os.path.join(WD, 'genes', 'all_genes_nonredundant')
    conda:
        "envs/filter_gtdbtk.yaml" 
    resources:
        cores = 40,
        memory = 188, 
        runtime = '1:00:00:00'
    shell:
        "bwa-mem2 index -p {output.index} {input.fasta}; samtools faidx {input.fasta}; touch {output.index}"

# Readmapping
rule bwa_readmap:
    input:
        index = os.path.join(WD, 'genes', 'all_genes_nonredundant'),
        fasta = os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta'),
        fastq1 = lambda wildcards: sample_dict[wildcards.sample][0],
        fastq2 = lambda wildcards: sample_dict[wildcards.sample][1]
    params:
        sample = SAMPLES
    output:
        bam = os.path.join(WD, 'mapped_reads', 'bams', 'gene_counts_{sample}.bam')
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores = 40,
        memory = 188,
        runtime = '1:00:00:00'
    shell:
        "bwa-mem2 mem -t {resources.cores} {input.index} {input.fastq1} {input.fastq2} | samtools view -T {input.fasta} -F 3584  -b --threads {resources.cores} | samtools sort --threads {resources.cores} > {output.bam}; samtools index {output.bam}"
