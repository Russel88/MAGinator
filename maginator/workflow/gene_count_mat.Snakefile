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

print(param_dict['map_filter'])
if param_dict['map_filter'] == 'pablo':
    rule filter_bamfile:
        input:
            os.path.join(WD,'mapped_reads', 'bams','gene_counts_{sample}.bam'),
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
            min_cov = param_dict['min_cov'],
            benchmark = param_dict['benchmark'],
        script:
            "scripts/extract_map.py"

else:
    rule filter_bamfile_shiraz:
        input:
            os.path.join(WD,'mapped_reads', 'bams','gene_counts_{sample}.bam'),
        output:
            os.path.join(WD,'mapped_reads', 'bams','gene_counts_{sample}_filtered.bam')
        conda:
            "envs/filter_geneclusters.yaml"
        resources:
            cores = 1,
            mem_gb = 20,
            runtime = 7200 #2h in s
        shell:
            """
            msamtools filter -b -l 80 -p 95 -z 80 {input} > {output}
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
rule filter_coverage:
    input:
        cov_file = os.path.join(WD,'mapped_reads','coverage','{sample}.coverage')         
    output:
        out_file = os.path.join(WD,'mapped_reads','counts','{sample}.filtered.counts'),
        names_file = os.path.join(WD,'mapped_reads','names','{sample}.gene_names')
    conda:
        "envs/filter_geneclusters.yaml"
    params:
        min_reads = param_dict['min_cov'],
        min_map = param_dict['min_map'],
        benchmark = param_dict['benchmark'],
        map_filter = param_dict['map_filter']
    resources:
        cores = 1,
        mem_gb = 20,  # Calculate the total size of input files in GB
        runtime =  7200 #2h in s
    script: 
        "scripts/coverage_filter.py"

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
    shell:
        "mv {input} {output}"

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
        mem_gb = 188,
        runtime = '3600' #1h in s
    shell:
        "paste {input.gene_names} {input.readcounts} | cat {input.header} - > {output}; sed -i '$d' {output}"