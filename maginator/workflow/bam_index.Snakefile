import os
import logging

WD = config['wd']
BAM_DIR = os.path.join(WD, config['bam_dir'])
SAMPLES = set(glob_wildcards(os.path.join(BAM_DIR, '{sample}.bam')).sample)

logging.info('Indexing BAM files')

rule all:
    input:
        expand(os.path.join(BAM_DIR, "{sample}.bam.bai"), sample=SAMPLES)

rule index:
    input:
        bam=os.path.join(BAM_DIR, "{sample}.bam")
    output:
        bai=os.path.join(BAM_DIR, "{sample}.bam.bai")
    shell:
        """
        samtools index {input.bam} {output.bai}
        """
