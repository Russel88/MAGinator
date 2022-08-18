import os
import re

WD = config['wd']
F_DIR = os.path.join(WD, config['fasta_dir'])
P_DIR = os.path.join(WD, config['prodigal_dir'])

CLUSTERS = set(glob_wildcards(os.path.join(F_DIR, '{cluster}')).cluster)
CLUSTERS = {x for x in CLUSTERS if x.isdigit()}

wildcard_constraints:
    cluster="\d+"

rule all:
    input:
        expand(os.path.join(P_DIR, 'nucl', '{cluster}.fnn'), cluster=CLUSTERS)

rule prodigal:
    input:
        os.path.join(F_DIR, '{cluster}/')
    output:
        faa = os.path.join(P_DIR, 'prot', '{cluster}.faa'),
        fnn = os.path.join(P_DIR, 'nucl', '{cluster}.fnn')
    conda:
        "envs/filter_prodigal.yaml"
    resources:
        cores=1,
        memory=4,
        runtime='02:00:00'
    shell:
        '''
        prodigal -i <(cat {input}/*.fa) -p meta -a {output.faa} -d {output.fnn} > /dev/null
        '''
