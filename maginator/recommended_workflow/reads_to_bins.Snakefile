#!python

#!/usr/bin/env python3
import os
import math
import Bio

#Get IDs / wildcards for the clusters and the samples

SAMPLES = set()
sample_dict = {}

with open(config['reads']) as f:
    for row in f:
        row = row.rstrip().split(',')
        SAMPLES.add(row[0])
        sample_dict[row[0]] = [row[1], row[2]]


rule all:
    input:
        expand("preprocessed/clean_{sample}_2.fastq", sample=SAMPLES),
        "vamb/clusters.tsv",
        "vamb/checkm.results",
        "preprocessed/read_stats.png"

# Unzipping the reads - outcomment if your data is not zipped
rule unzip:
    input:
        fastq1 = lambda wildcards: sample_dict[wildcards.sample][0],
        fastq2 = lambda wildcards: sample_dict[wildcards.sample][1]
    output:
        unzip1 = temp("data/{sample}_1.fastq"),
        unzip2 = temp("data/{sample}_2.fastq")
    params:
       sample = SAMPLES
    resources:
        cores = 1,
        memory = 45,
        runtime = '1:00:00'
    shell:
        "zcat {input.fastq1} > {output.unzip1}; zcat {input.fastq2} > {output.unzip2}"

# Removing adapters
rule adapter_removal:
    input:
        fastq1 = "data/{sample}_1.fastq",
        fastq2 = "data/{sample}_2.fastq"
    output:
        no_adapt1 = temp("preprocessed/no_adapt_{sample}_1.fastq"),
        no_adapt2 = temp("preprocessed/no_adapt_{sample}_2.fastq"),
        counts = temp("preprocessed/stats/counts_{sample}.tsv")
    params:
        sample = SAMPLES
    resources:
        cores = 1,
        memory = 62,
        runtime = '15:00:00'
    log:
        "log/read_qc/adapter_removal/{sample}.adapter.log"
    conda:
        "envs/preprocess.yaml"
    shell:
        "bbduk.sh in1={input.fastq1} in2={input.fastq2} out1={output.no_adapt1} out2={output.no_adapt2} ktrim=r k=23 mink=11 hdist=1 hdist2=0 forcetrimleft=3 tpe tbo 2>{log};"
        "mkdir -p $(dirname {output.counts}); echo $(cat {input.fastq1}|wc -l)/4|bc >> {output.counts};"
        "echo $(cat {output.no_adapt1}|wc -l)/4|bc >> {output.counts};"

if "read_length_cutoff" not in config:
    config["read_length_cutoff"] = 100

# Removing low quality reads
rule quality_removal:
    input:
        no_adapt1 = "preprocessed/no_adapt_{sample}_1.fastq",
        no_adapt2 = "preprocessed/no_adapt_{sample}_2.fastq",
        counts = "preprocessed/stats/counts_{sample}.tsv"
    output:
        qc_reads1 = temp("preprocessed/quality_{sample}_1.fastq"),
        qc_reads2 = temp("preprocessed/quality_{sample}_2.fastq"),
        qc_counts = temp("preprocessed/stats/qc_counts_{sample}.tsv")
    params:
        length = config["read_length_cutoff"]
    resources:
        cores = 1,
        memory = 45,
        runtime = '5:00:00'
    log:
        "log/read_qc/quality_filtering/{sample}.quality.log"
    conda:
        "envs/preprocess.yaml"
    shell:
        "sickle pe --qual-type sanger -f {input.no_adapt1} -r {input.no_adapt2} -o {output.qc_reads1} -p {output.qc_reads2} -s /dev/null -l {params.length} 2>{log};"
        "cp {input.counts} {output.qc_counts}; echo $(cat {output.qc_reads1}|wc -l)/4|bc >> {output.qc_counts};"

# Import the human genome
rule hg19:
#    input:
#        dummy = "preprocessed/empty.txt"
    output:
        hg19 = "preprocessed/hg19.fasta"
    resources:
        cores = 1,
        memory = 45,
        runtime = '5:00:00'
    conda:
        "envs/preprocess.yaml"
    script:
        "envs/import_hg_19.R"

# Indexing the human genome
rule human_index:
    input: 
        hg19 = ancient("preprocessed/hg19.fasta")
    output:
        ref_dir = directory("preprocessed/ref"),
        dummy = "preprocessed/empty.txt"
    resources:
        cores = 20,
        memory = 179,
        runtime = '20:00:00'
    log:
        "log/read_qc/human_removal/index_human.log"
    conda:
        "envs/preprocess.yaml"
    shell:
        "bbmap.sh ref={input.hg19} 2>{log}; mv ref/ {output.ref_dir}; touch {output.dummy}"

# Removing human contamination
rule human_removal:
    input:
        qc_reads1 = "preprocessed/quality_{sample}_1.fastq",
        qc_reads2 = "preprocessed/quality_{sample}_2.fastq",
        dummy = ancient("preprocessed/empty.txt"),
        qc_counts = "preprocessed/stats/qc_counts_{sample}.tsv"
    output:
        clean_reads1 = "preprocessed/clean_{sample}_1.fastq",
        clean_reads2 = "preprocessed/clean_{sample}_2.fastq",
        human_counts = temp("preprocessed/stats/human_counts_{sample}.tsv")
    resources:
        cores = 5,
        memory = 180,
        runtime = '25:00:00:00'
    log:
        "log/read_qc/human_removal/{sample}.human.log"
    conda:
        "envs/preprocess.yaml"
    shell:
        "bbmap.sh qin=33 threads={resources.cores} in={input.qc_reads1} in2={input.qc_reads2} outu1={output.clean_reads1} outu2={output.clean_reads2} path=preprocessed -Xmx{resources.memory}g 2>{log};"
        "cp {input.qc_counts} {output.human_counts}; echo $(cat {output.clean_reads1}|wc -l)/4|bc >> {output.human_counts};"

# Combining the read statistics
rule stat_combine:
    input:
        expand("preprocessed/stats/human_counts_{sample}.tsv", sample=SAMPLES)
    output:
        "preprocessed/stats/read_stats.tsv"
    resources:
        cores = 1,
        memory = 20,
        runtime = '1:00:00'
    log:
        "log/read_qc/read_counts_combine.log"
    shell:
        "paste {input} > {output} 2>{log};" 

# creating reads stats figure
rule stats_fig:
    input:
        stats = "preprocessed/stats/read_stats.tsv"
    output:
        fig = "preprocessed/read_stats.png"
    resources:
        cores = 1,
        memory = 20,
        runtime = '1:00:00'
    run:
        import numpy as np
        from matplotlib import pyplot as plt

        a = np.loadtxt(input.stats, delimiter='\t')
        A = a[::, a[3,].argsort()[::-1]] # sorting the np array according to amount of reads in the sample
        reads_left = A[3]
        human_reads = A[2]-A[3]
        low_qc_reads = A[1]-A[2]
        adapter_reads = A[0]-A[1]

        ind = np.arange(len(A[0])) + 0.75

        fig, ax = plt.subplots(1,1)
        p0 = ax.bar(ind, reads_left, color = 'lightblue')
        p1 = ax.bar(ind, adapter_reads, bottom = reads_left, color = 'blue')
        p2 = ax.bar(ind, low_qc_reads, bottom = reads_left + adapter_reads, color = 'grey')
        p3 = ax.bar(ind, human_reads, bottom = reads_left + adapter_reads + low_qc_reads, color = 'r')

        plt.xticks([], [])
        ax.set_ylabel('Reads')
        ax.set_xlabel('Samples')
        fig.legend( (p0[0], p1[0], p2[0], p3[0]), ('Reads', 'Adapter reads', 'Low quality reads', 'Human reads') )
        fig.savefig(output.fig)


# Metagenomic assembly of each sample
rule assembly:
    input:
        R1 = "preprocessed/clean_{sample}_1.fastq",
        R2 = "preprocessed/clean_{sample}_2.fastq"
    output:
        "assembly/{sample}/contigs.fasta"
    resources:
        cores = 40,
        memory = 185,
        runtime = '20:00:00:00'
    params:
        statusfile = "assembly/{sample}/pipeline_state/stage_0_before_start"
    log:
        "log/assembly/metaspades.{sample}.log"
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        # find path for output .paths file - will indicate whether spades has been run before and can be rerun.  
        if [ -f {params.statusfile} ]; then (metaspades.py -o $(dirname {output}) --continue 2>{log})
        else (metaspades.py -t {resources.cores} -m {resources.memory} -o $(dirname {output}) -1 {input.R1} -2 {input.R2} -k 21,33,55,77,99,127 2>{log} )
        fi
        """

# Add the sampleID to the contigID  
rule contig_rename:
    input:
        "assembly/{sample}/contigs.fasta"
    output:
        "assembly/{sample}/renamed_contigs.fasta"
    resources:
        cores = 1,
        memory = 45,
        runtime = '5:00:00'
    log:
        "log/assembly/{sample}.id.log"
    shell:
        "cat {input} | sed 's/>/>{wildcards.sample}@/g' > {output} 2>{log}"

# Combining all the contigs
rule combine_contigs:
    input:
        expand("assembly/{sample}/renamed_contigs.fasta", sample=SAMPLES)
    output:
        "assembly/all_assemblies_all_lengths.fasta"
    resources:
        cores = 1,
        memory = 45,
        runtime = '5:00:00'
    log:
        "log/assembly/combine_contigs.log"
    shell:
        "cat {input} >> {output} 2>{log}"


# sorting away contigs less than contig_length_cutoff (default: 1500bp)
if "contig_length_cutoff" not in config:
    config["contig_length_cutoff"] = 1500

rule contig_filtering:
    input:
        all_lengths = "assembly/all_assemblies_all_lengths.fasta"
    output:
        filtered_lengths = "assembly/all_assemblies.fasta"
    params:
      size = config["contig_length_cutoff"]
    resources:
        cores = 4,
        memory = 180,
        runtime = '15:00:00'
    run:
        from Bio import SeqIO
        long_sequences = []
        with open(input.all_lengths, "r") as handle:

            for record in SeqIO.parse(handle, "fasta"):
                if len(record) >= int(params.size):
                    long_sequences.append(record)

        SeqIO.write(long_sequences, output.filtered_lengths, "fasta")


# index the repressentative nucleotide file
rule index:
    input:
       "assembly/all_assemblies.fasta"
    output:
       done=touch("assembly/all_assemblies"),
       index="assembly/all_assemblies.fasta.bwt"
    resources:
        memory = 1450,
        runtime = '2:00:00:00',
        cores = 40
    log:
        "log/mapping/index.log"
    conda:
        "envs/preprocess.yaml"
    shell:
        "bwa-mem2 index {input}; touch {output.index} 2>{log}"

# map the repressentative sequences back to the original files
rule map:
    input:
        index="assembly/all_assemblies",
        gene_cat="assembly/all_assemblies.fasta",
        R1 = lambda wildcards: sample_dict[wildcards.sample][0],
        R2 = lambda wildcards: sample_dict[wildcards.sample][1]
    output:
        "mapped/{sample}.bam"
    params:
        sample=SAMPLES
    resources:
        cores = 14,
        memory = 100,
        runtime = '1:00:00:00'
    log:
        "log/mapping/map_{sample}.log"
    conda:
        "envs/preprocess.yaml"
    shell:
        "module unload samtools; module load samtools/1.10; bwa-mem2 mem -t 14 {input.gene_cat} {input.R1} {input.R2} | samtools view -Sb -F4 - > {output} 2>{log}"



#####################################################################################################################################################################################################################################
########################################################## The rest of the workflow has been copied from VAMB (https://github.com/RasmussenLab/vamb/blob/master/workflow/) ##########################################################
#####################################################################################################################################################################################################################################

rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        temp("mapped/{sample}.sort.bam")
    resources:
        cores = 1,
        memory = 15,
        runtime = '10:00:00:00'
    params:
        prefix="mapped/tmp.{sample}"
    log:
        "log/map/{sample}.sort.log"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort {input} -T {params.prefix} --threads {resources.cores} -m 3G -o {output} 2>{log}"

rule jgi:
    input:
        bam = "mapped/{sample}.sort.bam"
    output:
        jgi = temp("jgi/{sample}.raw.jgi")
    resources:
        cores = 1,
        memory = 10,
        runtime = '10:00:00:00'
    log:
        "log/jgi/{sample}.jgi"
    conda:
        "envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths --noIntraDepthVariance --outputDepth {output} {input} 2>{log}"

rule cut_column1to3: 
    input:
        "jgi/%s.raw.jgi" % SAMPLES.pop() 
    output:
        "jgi/jgi.column1to3"
    resources:
        cores = 1,
        memory = 1,
        runtime = '1:00:00:00'
    log:
        "log/jgi/column1to3"
    shell:
        "cut -f1-3 {input} > {output} 2>{log}"

rule cut_column4to5:
    input:
        "jgi/{sample}.raw.jgi"
    output:
        "jgi/{sample}.cut.jgi"
    resources:
        cores = 1,
        memory = 1,
        runtime = '1:00:00:00'
    log:
        "log/jgi/{sample}.cut.log"
    shell:
        "cut -f1-3 --complement {input} > {output} 2>{log}"

rule paste_abundances:
    input:
        column1to3="jgi/jgi.column1to3",
        data=expand("jgi/{sample}.cut.jgi", sample=SAMPLES)
    output:
        "jgi_matrix/jgi.abundance.dat"
    resources:
        cores = 1,
        memory = 1,
        runtime = '1:00:00:00'
    log:
        "log/jgi/paste_abundances.log"
    shell:
        "paste {input.column1to3} {input.data} > {output} 2>{log}"


rule vamb:
    input:
        jgi = "jgi_matrix/jgi.abundance.dat",
        contigs = "assembly/all_assemblies.fasta"
    output:
        "vamb/clusters.tsv",
        "vamb/latent.npz",
        "vamb/lengths.npz",
        "vamb/log.txt",
        "vamb/model.pt",
        "vamb/mask.npz",
        "vamb/tnf.npz"
    resources:
        cores = 20,
        memory = 180,
        runtime = '1:00:00:00'
    log:
        "log/vamb/vamb.log"
    conda:
        "envs/vamb.yaml"
    shell:
        "rm -rf vamb;"
        "vamb --outdir vamb --fasta {input.contigs} --jgi {input.jgi} -o @ -m 2000 --minfasta 2500 2>{log}"

rule checkm:
    input:
        "vamb/clusters.tsv"
    output:
        "vamb/checkm.results"
    resources:
        cores = 10,
        memory = 25,
        runtime = '15:00:00:00'
    params:
        bins = "vamb/bins",
        outdir = "vamb/checkm.outdir"
    log:
        "log/vamb/checkm.log"
    conda:
        "envs/checkm.yaml"
    shell:
        "checkm lineage_wf -f {output} -t {resources.cores} -x fna {params.bins} {params.outdir} 2>{log}"
