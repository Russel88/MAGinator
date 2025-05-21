import os
import re
import subprocess
import math

WD = config['wd']
PARAMS = config['params']
READS = config['reads']
with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}


rule all:
    input:
        os.path.join(WD, 'gtdbtk', 'gtdbtk.done')


rule gtdbtk:
    input:
        directory(os.path.join(WD,'combined_clusters'))
    output:
        directory(os.path.join(WD, 'gtdbtk')),
        os.path.join(WD, 'gtdbtk', 'gtdbtk.done')
    conda:
        "envs/filter_gtdbtk.yaml"
    params:
        gtdbtk=param_dict['gtdb_db']
    resources:
        cores=10,
        mem_gb=180,
        runtime=24*60*60*60
    shell:
        '''
        export GTDBTK_DATA_PATH={params.gtdbtk};
        gtdbtk classify_wf --genome_dir {input} --out_dir {output[0]} --cpus {resources.cores} --extension fa --keep_intermediates --skip_ani_screen
        if [[ -f {output[0]}/classify/gtdbtk.bac120.summary.tsv || -f {output[0]}/classify/gtdbtk.ar53.summary.tsv ]]; then
            touch {output[1]}
        fi
        '''
