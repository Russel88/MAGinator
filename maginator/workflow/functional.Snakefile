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
        os.path.join(WD, 'abundance', 'functional_profile.RData')

rule annotation:
    input:
        os.path.join(WD, 'genes', 'all_genes.faa')
    output:
        os.path.join(WD, 'functional', 'all_genes_annotation.tsv')
    conda:
        "envs/functional.yaml"
    resources:
        cores = 40,
        mem_gb = 40,
        runtime = 72000 #20h in s
    script:
        "emapper.py -i {input} -o {output} -m diamond --cpu {resources.cores}"

# Merging the annotation with the gene count matrix into an RData object
rule functional_profiling:
    input:
        os.path.join(WD,  'functional', 'all_genes_annotation.tsv')
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')
    output:
        os.path.join(WD, 'abundance', 'functional_profile.RData')
    conda:
        "envs/functional.yaml"
    resources:
        cores = 1,
        mem_gb = 40,
        runtime = 7200 #2h in s
    script:
        "scripts/functional_profiling.R"