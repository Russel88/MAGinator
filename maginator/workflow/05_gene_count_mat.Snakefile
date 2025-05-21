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


rule bam_name_sorting:
    input:
        bam = os.path.join(WD, 'mapped_reads', 'bams', 'gene_counts_{sample}.bam')
    output:
        bam = os.path.join(WD, 'mapped_reads', 'bams', 'gene_counts_{sample}.sorted.bam')
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores = 40,
        mem_gb = 188,
        runtime = 86400 #1d in s
    shell:
        "samtools sort --threads {resources.cores} -o {output.bam} {input.bam}"


rule filter_bamfile:
    input:
        os.path.join(WD,'mapped_reads', 'bams','gene_counts_{sample}.sorted.bam'),
    output:
        os.path.join(WD,'mapped_reads', 'bams','gene_counts_{sample}_filtered.bam')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 7200 #2h in s
    params: 
        min_map = param_dict['min_map'],
        min_iden = param_dict['min_identity'],
        min_len = param_dict['min_length'],
    shell:
        """
        msamtools filter -b -l {params.min_len} -p {params.min_iden} -z {params.min_map} {input} > {output}
        """  

# Use samtools to count the genes for each sample
rule gene_coverage:
    input:
        os.path.join(WD,'mapped_reads', 'bams','gene_counts_{sample}_filtered.bam'),
    output:
        os.path.join(WD,'mapped_reads','coverage','{sample}.coverage')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 7200 #2h in s
    shell:
        "samtools coverage {input} > {output}"

#Modifying the gene count matrix so genes that do not reach a threshold of mapped counts are set to 0   (Find a way to include the gene name as output)
rule extract_coverage:
    input:
        cov_file = os.path.join(WD,'mapped_reads','coverage','{sample}.coverage')         
    output:
        out_file = os.path.join(WD,'mapped_reads','counts','{sample}.filtered.counts'),
        names_file = os.path.join(WD,'mapped_reads','names','{sample}.gene_names')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 20,  # Calculate the total size of input files in GB
        runtime =  7200 #2h in s
    script: 
        "scripts/extract_coverage.py"

# Get only the names from one file, and find a way to safely delete the rest (they are the same)
rule clear_gene_names:
    input:
        os.path.join(WD, 'mapped_reads', 'names', '{sample}.gene_names')
    output:
        os.path.join(WD, 'mapped_reads', 'gene_names_{sample}')
    wildcard_constraints:
        sample=SAMPLES[0]
    resources:
        cores = 1,
        mem_gb = 5,
        runtime = '3600' #1h in s
    params:
        names_dir = os.path.join(WD, 'mapped_reads', 'names')    
    shell:
        "mv {input} {output}; rm -r {params.names_dir}"

# Create the header for the gene count matrix
rule create_header:
    input: 
        expand(os.path.join(WD, 'mapped_reads', 'counts', '{sample}.filtered.counts'), sample=SAMPLES)
    output:
        header = os.path.join(WD, 'mapped_reads', 'header.txt')
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = '3600' #1h in s
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
        readcounts = expand(os.path.join(WD, 'mapped_reads', 'counts', '{sample}.filtered.counts'), sample=SAMPLES),
        gene_names = expand(os.path.join(WD, 'mapped_reads', 'gene_names_{sample}'), sample=SAMPLES[0]),
        header = os.path.join(WD, 'mapped_reads', 'header.txt')
    output:
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 40,
        runtime = '216000' #1h
    shell:
        "paste {input.gene_names} {input.readcounts} | cat {input.header} - > {output}; sed -i '$d' {output}"
