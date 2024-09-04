import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        read_counts = os.path.join(WD, 'benchmark', 'mapped_reads.tsv'),
        gene_counts = os.path.join(WD, 'benchmark', 'genes_with_reads.tsv')

rule benchmarking:
    input:
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')
    output:
        read_counts = os.path.join(WD, 'benchmark', 'mapped_reads.tsv'),
        gene_counts = os.path.join(WD, 'benchmark', 'genes_with_reads.tsv')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 40,
        runtime = '3600' #1h in s
    params:
        min_map = param_dict['min_map'],
        min_iden = param_dict['min_identity'],
        min_len = param_dict['min_length'],
    script:
        "scripts/benchmark.py"