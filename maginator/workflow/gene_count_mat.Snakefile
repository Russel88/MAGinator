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
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')


# Use samtools to count the genes for each sample
rule count_genes:
    input:
        os.path.join(WD, 'mapped_reads', 'bams', 'gene_counts_{sample}.bam')
    output:
        readcounts = os.path.join(WD, 'mapped_reads', 'counts', '{sample}.counts')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        memory = 20,
        runtime = '2:00:00'
    shell:
        "samtools idxstats {input} | cut -f3 > {output.readcounts}"

# Create the gene index
rule gene_names:
    input:
        os.path.join(WD, 'mapped_reads', 'bams', 'gene_counts_{sample}.bam')
    output:
        os.path.join(WD, 'mapped_reads', 'gene_names_{sample}')
    wildcard_constraints:
        sample=SAMPLES[0]
    conda:
       	"envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        memory = 20,
        runtime = '1:00:00'
    shell:
        "samtools idxstats {input} | cut -f1 > {output}"


# Create the header for the gene count matrix
rule create_header:
    input: 
        expand(os.path.join(WD, 'mapped_reads', 'counts', '{sample}.counts'), sample=SAMPLES)
    output:
        header = os.path.join(WD, 'mapped_reads', 'header.txt')
    resources:
        cores = 1,
        memory = 188,
        runtime = '1:00:00'
    run:
        header = "Gene"
        for f in input:
            sample_name = f.split("/")[-1].split(".")[0]
            header = header + "\t" + sample_name
        header = header + "\n"
        with open(output.header, "w") as out:
            out.write(header)


#Combining the gene names with the counts of all samples
rule gene_count_matrix:
    input:
        readcounts = expand(os.path.join(WD, 'mapped_reads', 'counts', '{sample}.counts'), sample=SAMPLES),
        gene_names = expand(os.path.join(WD, 'mapped_reads', 'gene_names_{sample}'), sample=SAMPLES[0]),
        header = os.path.join(WD, 'mapped_reads', 'header.txt')
    output:
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')
    conda:
       	"envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        memory = 188,
        runtime = '1:00:00'
    shell:
        "paste {input.gene_names} {input.readcounts} | cat {input.header} - > {output}; sed -i '$d' {output}"
