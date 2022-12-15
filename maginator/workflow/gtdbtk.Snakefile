import os
import re
import subprocess
import math

WD = config['wd']
F_DIR = os.path.join(WD, 'clusters')
PARAMS = config['params']
READS = config['reads']
with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

CLUSTERS = set(glob_wildcards(os.path.join(F_DIR, '{cluster}')).cluster)
CLUSTERS = {x for x in CLUSTERS if x.isdigit()}

# Adaptive resource usage
## Get lines in read file
out = subprocess.Popen(['wc', '-l', READS], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
n_samples = int(out.partition(b' ')[0])
## time is 1 hour per 30 samples
tim = str(math.ceil(n_samples/30))+':00:00'

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
        runtime=tim
    shell:
        '''
        export GTDBTK_DATA_PATH={params.gtdbtk};
        gtdbtk classify_wf --genome_dir {input} --out_dir {output} --cpus {resources.cores} --extension fa --keep_intermediates || true
        '''
