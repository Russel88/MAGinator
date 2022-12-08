import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

VAMB = config['vamb']

rule all:
    input:
        os.path.join(WD, 'tabs', 'gene_cluster_tax_scope.tab'),
        os.path.join(WD, 'tabs', 'gene_cluster_bins.tab')

rule gene_tax:
    input:
        tax=os.path.join(WD, 'tabs', 'tax_matrix.tsv'),
        mgs=os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        cluster=os.path.join(WD, 'genes', 'all_genes_cluster.tsv'),
        vamb=VAMB
    output:
        os.path.join(WD, 'tabs', 'gene_cluster_tax_scope.tab'),
        os.path.join(WD, 'tabs', 'gene_cluster_bins.tab')
    params:
        cutoff = param_dict['tax_scope_threshold']
    conda:
        "envs/phylo.yaml"
    resources:
        cores = 1,
        memory = 10,
        runtime = '02:00:00'
    script:
        "scripts/gene_cluster2tax.py"

