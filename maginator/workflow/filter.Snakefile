import os

WD = config['wd']
CONTIGS = config['contigs']
VAMB = config['vamb']
F_DIR = os.path.join(WD, config['fasta_dir'])

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
        binsize=config['binsize']
    conda:
        "envs/filter_prodigal.yaml"
    resources:
        cores=1,
        memory=32,
        runtime='02:00:00'
    script:
        "../filter.py"

