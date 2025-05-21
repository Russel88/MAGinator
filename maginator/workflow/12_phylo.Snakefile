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

# Fasttree
if param_dict['phylo'] == 'fasttree':
    rule tree:
        input:
            os.path.join(WD, 'phylo', 'cluster_alignments', '{cluster}.fna')
        output:
            os.path.join(WD, 'phylo', 'trees', '{cluster}.nwk')
        conda:
            "envs/phylo.yaml"
        resources:
            cores=1,
            mem_gb=32,
            runtime=43200 #12h in s
        shell:
            """
            fasttree -nt {input} > {output}
            """

# IQtree
if param_dict['phylo'] == 'iqtree':
    rule tree:
        input:
            fna=os.path.join(WD, 'phylo', 'cluster_alignments', '{cluster}.fna'),
            part=os.path.join(WD, 'phylo', 'cluster_alignments', '{cluster}.part')
        output:
            os.path.join(WD, 'phylo', 'trees', '{cluster}.nwk')
        conda:
            "envs/phylo.yaml"
        params:
            prefix=os.path.join(WD, 'phylo', 'trees', '{cluster}'),
        resources:
            cores=40,
            mem_gb=180,
            runtime=172800 #2d in s
        shell:
            """
            iqtree -T {resources.cores} -s {input.fna} -p {input.part} -o Outgroup --prefix {params.prefix} || true && if [ -f {params.prefix}.treefile ]; then mv {params.prefix}.treefile {output}; else touch {output}; fi
            """
