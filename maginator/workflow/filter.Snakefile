import os
import subprocess
import math
import distutils.util

WD = config['wd']
CONTIGS = config['contigs']
VAMB = config['vamb']
F_DIR = os.path.join(WD, 'clusters')
PARAMS = config['params']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

# Adaptive resource usage
## Get lines in vamb file
out = subprocess.Popen(['wc', '-l', VAMB], stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
n_contigs = int(out.partition(b' ')[0])
## mem_gb is 30gb per million contigs
mem = math.ceil(n_contigs/1000000)*30
if mem > int(param_dict['max_mem']):
    mem = int(param_dict['max_mem'])
## time is 1 hour per million
tim = str(max(1800, min(math.ceil(n_contigs/100000)*60*60, 72000))) # runtime in seconds, max 20h, min 30min

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
        mem_gb=mem,
        runtime=tim
    script:
        "scripts/filter.py"

