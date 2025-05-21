import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

with open(param_dict['reads']) as f:
    SAMPLES = ([row.split(',')[0] for row in f])

rule all:
    input:
        expand(os.path.join(WD, 'phylo', 'pileup', '{sample}.mp'), sample=SAMPLES)

rule subset:
    input:
        bed=os.path.join(WD, 'phylo', 'intermediate', 'cluster_genes.bed'),
        bam=os.path.join(WD, 'mapped_reads', 'bams', 'gene_counts_{sample}.bam')
    output:
        os.path.join(WD, 'mapped_reads', 'bams', 'subset', '{sample}.bam')
    conda:
        "envs/phylo.yaml"
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 43200 #12h in s
    shell:
        "samtools view -b -L {input.bed} {input.bam} | samtools sort -o {output}"

rule pileup:
    input:
        fna=os.path.join(WD, 'phylo', 'intermediate', 'cluster_genes.fna'),
        bam=os.path.join(WD, 'mapped_reads', 'bams', 'subset', '{sample}.bam')
    output:
        os.path.join(WD, 'phylo', 'pileup', '{sample}.mp')
    conda:
        "envs/phylo.yaml"
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 86400 #1d in s
    shell:
        "samtools mpileup -A -x -f {input.fna} {input.bam} -o {output}"

