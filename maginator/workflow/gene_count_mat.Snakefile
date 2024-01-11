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

rule filter_bamfile:
    input:
        os.path.join(WD,'mapped_reads', 'bams','gene_counts_{samples}.bam'),
    output:
        os.path.join(WD,'mapped_reads', 'bams','gene_counts_{samples}_filtered.bam')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 7200 #2h in s
    params:
        min_map = param_dict['min_map']
    script:
        "scripts/extract_map.py"  

# Use samtools to count the genes for each sample
rule gene_coverage:
    input:
        os.path.join(WD,'mapped_reads', 'bams','gene_counts_{samples}_filtered.bam'),
    output:
        os.path.join(WD,'mapped_reads','coverage','{samples}.coverage')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 7200 #2h in s
    shell:
        "samtools coverage {input} > {output}"

#Modifying the gene count matrix so genes that do not reach a threshold of mapped counts are set to 0  
rule create_gene_count_matrix:
    input:
        cov_file = os.path.join(WD,'mapped_reads','coverage','{samples}.coverage')         
    output:
        out_file = os.path.join(WD,'mapped_reads','coverage','{samples}_filtered.coverage')
        gene_names = os.path.join(WD, 'mapped_reads', 'gene_names')
    conda:
        "envs/filter_geneclusters.yaml"
    params:
        min_reads = param_dict['min_cov'],
            
    resources:
        cores = 1,
        mem_gb = 20,  # Calculate the total size of input files in GB
        runtime =  7200 #2h in s
    script: 
        "scripts/coverage_filter.py"

# Create the header for the gene count matrix
rule create_header:
    input: 
        expand(os.path.join(WD, 'mapped_reads', 'counts', '{sample}_filtered.counts'), sample=SAMPLES)
    output:
        header = os.path.join(WD, 'mapped_reads', 'header.txt')
    resources:
        cores = 1,
        mem_gb = 188,
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
        readcounts = expand(os.path.join(WD, 'mapped_reads', 'counts', '{sample}_filtered.counts'), sample=SAMPLES),
        gene_names = os.path.join(WD, 'mapped_reads', 'gene_names'),
        header = os.path.join(WD, 'mapped_reads', 'header.txt')
    output:
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')
    conda:
        "envs/filter_geneclusters.yaml"
    resources:
        cores = 1,
        mem_gb = 188,
        runtime = '3600' #1h in s
    shell:
        "paste {input.gene_names} {input.readcounts} | cat {input.header} - > {output}; sed -i '$d' {output}"