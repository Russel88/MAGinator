import os
import re

WD = config['wd']
PARAMS = config['params']
with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

CLUSTERS = set(glob_wildcards(os.path.join(WD, 'phylo', 'cluster_alignments', '{cluster}.fna')).cluster)

rule all:
    input:
        expand(os.path.join(WD, 'phylo', 'trees', '{cluster}.nwk'), cluster=CLUSTERS),

rule fasttree:
    input:
        os.path.join(WD, 'phylo', 'cluster_alignments', '{cluster}.fna')
    output:
        os.path.join(WD, 'phylo', 'trees', '{cluster}.nwk')
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        memory=32,
        runtime='12:00:00'
    shell:
        """
        fasttree -nt {input} > {output}
        """

