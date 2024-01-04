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
## time is 1 hour per 30 samples or at least 2h
tim = (max(math.ceil(n_samples/30)*60, 2*60)) # time in h

wildcard_constraints:
    cluster="\d+"

rule all:
    input:
        expand(os.path.join(WD, 'gtdbtk', '{cluster}/gtdbtk.done'), cluster=CLUSTERS)

rule gtdbtk:
    input:
        ancient(os.path.join(F_DIR, '{cluster}/'))
    output:
        directory(os.path.join(WD, 'gtdbtk', '{cluster}')),
        os.path.join(WD, 'gtdbtk', '{cluster}/gtdbtk.done')
    conda:
        "envs/filter_gtdbtk.yaml"
    params:
        gtdbtk=param_dict['gtdb_db']
    resources:
        cores=8,
        mem_gb=180,
        runtime=str(tim*60*60)
    shell:
        '''
        export GTDBTK_DATA_PATH={params.gtdbtk};
        gtdbtk classify_wf --genome_dir {input} --out_dir {output[0]} --cpus {resources.cores} --extension fa --keep_intermediates --skip_ani_screen || true
        if [[ -f {output[0]}/classify/gtdbtk.bac120.summary.tsv || -f {output[0]}/classify/gtdbtk.ar53.summary.tsv ]]; then
            touch {output[1]}
        fi
        '''
