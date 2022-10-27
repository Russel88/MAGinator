import os

WD = config['wd']
CONTIGS = config['contigs']
VAMB = config['vamb']
F_DIR = os.path.join(WD, 'clusters')
PARAMS = config['params']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        F_DIR

rule bin_filter:
    input:
        VAMB,
        CONTIGS
    output:
        directory(F_DIR),
    params:
        binsize=param_dict['binsize']
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores=1,
        memory=32,
        runtime='02:00:00'
    script:
        "scripts/filter.py"

