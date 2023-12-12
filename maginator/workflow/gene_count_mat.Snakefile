import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}


with open(param_dict['reads']) as f:
    SAMPLES = ([row.split(',')[0] for row in f])

rule all:
    input:
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')

# Use samtools to count the genes for each sample
rule gene_coverage:
    input:
        os.path.join(WD,'mapped_reads', 'bams','gene_counts_{samples}.bam'),
    output:
        os.path.join(WD,'mapped_reads','coverage','{samples}.coverage')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 7200 #2h in s
    shell:
        "samtools coverage {input} > {output}"

#Modifying the gene count matrix so genes that do not reach a threshold of mapped counts are set to 0  
rule create_gene_count_matrix:
    input:
        cov_files = expand(os.path.join(WD,'mapped_reads','coverage','{samples}.coverage'),samples=SAMPLES)         
    output:
        gene_matrix=os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')
    conda:
        "envs/signature_genes.yaml"
    params:
        min_reads = param_dict['min_cov'],
    resources:
        cores = 1,
        mem_gb = 40,
        runtime = 7200 #2h in s
    script: 
        "scripts/gene_count_matrix.R"