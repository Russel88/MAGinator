import os
import re

WD = config['wd']
PARAMS = config['params']
with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

CLUSTERS = set(glob_wildcards(os.path.join(WD, 'phylo', 'fastas', '{cluster}', '{fastas}')).cluster)

rule all:
    input:
        expand(os.path.join(WD, 'phylo', 'alignments', '{cluster}'), cluster=CLUSTERS),
        os.path.join(WD, 'phylo', 'cluster_alignments'),

rule align:
    input:
        os.path.join(WD, 'phylo', 'fastas', '{cluster}')
    output:
        directory(os.path.join(WD, 'phylo', 'alignments', '{cluster}'))
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        mem_gb=4,
        runtime=36000 #10h in s
    shell:
        """
        mkdir -p {output}
        for gene in {input}/*.fna; do mafft --ep 0.123 $gene > {output}/$(basename $gene); done
        """

rule concat:
    input:
        param_dict['reads'],
        expand(os.path.join(WD, 'phylo', 'alignments', '{cluster}'), cluster=CLUSTERS)
    output:
        directory(os.path.join(WD, 'phylo', 'cluster_alignments'))
    conda:
        "envs/phylo.yaml"
    resources:
        cores=10,
        mem_gb=40,
        runtime=36000 #10h in s
    params:
        ali_dir=os.path.join(WD, 'phylo', 'alignments')
    script:
        "scripts/concat.py"


