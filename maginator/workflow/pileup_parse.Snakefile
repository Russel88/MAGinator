import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        os.path.join(WD, 'phylo', 'fastas'),
        os.path.join(WD, 'phylo', 'stats.tab'),
        os.path.join(WD, 'phylo', 'stats_genes.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'pileup_dumb.bin')

rule parse:
    input:
        pileup=os.path.join(WD, 'phylo', 'pileup'),
        ref=os.path.join(WD, 'phylo', 'intermediate', 'cluster_genes.fna'),
        root=os.path.join(WD, 'phylo', 'intermediate', 'markers_roots.fna'),
        sig=os.path.join(WD, 'tabs', 'signature_genes_cluster.tsv'),
        mark=os.path.join(WD, 'phylo', 'intermediate', 'clusterID_markerGene_geneCluster_rootGene.tsv'),
        tax=os.path.join(WD, 'tabs', 'tax_matrix.tsv')
    output:
        align=directory(os.path.join(WD, 'phylo', 'fastas')),
        stat=os.path.join(WD, 'phylo', 'stats.tab'),
        gene=os.path.join(WD, 'phylo', 'stats_genes.tab'),
        dump=os.path.join(WD, 'phylo', 'intermediate', 'pileup_dumb.bin')
    params:
        min_af = param_dict['min_af'], 
        min_depth = param_dict['min_depth'], 
        min_nonN = param_dict['min_nonN'], 
        min_marker_genes = param_dict['min_marker_genes'], 
        min_signature_genes = param_dict['min_signature_genes']
    conda:
        "envs/phylo.yaml"
    resources:
        cores = 40,
        memory = 180,
        runtime = '02:00:00:00'
    script:
        "scripts/mpileup.py"

