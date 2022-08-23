import os
import re

WD = config['wd']
F_DIR = os.path.join(WD, 'clusters')
PARAMS = config['params']
with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

CLUSTERS = set(glob_wildcards(os.path.join(F_DIR, '{cluster}')).cluster)
CLUSTERS = {x for x in CLUSTERS if x.isdigit()}

wildcard_constraints:
    cluster="\d+"

rule all:
    input:
        expand(os.path.join(WD, 'gtdbtk', '{cluster}'), cluster=CLUSTERS)

rule gtdbtk:
    input:
        os.path.join(F_DIR, '{cluster}/')
    output:
        directory(os.path.join(WD, 'gtdbtk', '{cluster}'))
    conda:
        "envs/filter_gtdbtk.yaml"
    params:
        gtdbtk=param_dict['gtdb_db']
    resources:
        cores=20,
        memory=90,
        runtime='02:00:00'
    shell:
        '''
        export GTDBTK_DATA_PATH={params.gtdbtk};
        gtdbtk classify_wf --genome_dir {input} --out_dir {output} --cpus {resources.cores} --extension fa --keep_intermediates || true
        '''
