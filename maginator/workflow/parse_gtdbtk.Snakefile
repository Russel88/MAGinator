import os

WD = config['wd']
PARAMS = config['params']
VAMB = config['vamb']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'genes', 'all_genes.fna'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_unique_bac_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_unique_ar_markers.tab')

rule parse_gtdbtk:
    input:
        os.path.join(WD, 'gtdbtk'),
        WD
    output:
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'genes', 'all_genes.fna'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_unique_bac_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_unique_ar_markers.tab')
    params:
        param_dict['annotation_prevalence'],
        VAMB
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores=1,
        memory=32,
        runtime='10:00:00'
    script:
        "../parse_gtdbtk.py"

