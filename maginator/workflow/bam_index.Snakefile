import os

WD = config['wd']
BAM_DIR = os.path.join(WD, config['bam_dir'])
SAMPLES = set(glob_wildcards(os.path.join(BAM_DIR, '{sample}.bam')).sample)

rule all:
    input:
        expand(os.path.join(BAM_DIR, "{sample}.bam.bai"), sample=SAMPLES)

rule index:
    input:
        bam=os.path.join(BAM_DIR, "{sample}.bam")
    output:
        bai=os.path.join(BAM_DIR, "{sample}.bam.bai")
    log:
        os.path.join(WD, "logs/bam_index/{sample}.log")
    shell:
        """
        samtools index {input.bam} {output.bai} > {log} 2>&1
        """
