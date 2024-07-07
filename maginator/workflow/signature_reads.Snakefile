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
        bam = os.path.join(WD, 'mapped_reads', 'bams', 'gene_counts_{sample}.name_sorted.bam')
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores = 40,
        mem_gb = 188,
        runtime = 86400 #1d in s
    shell:
        "samtools sort -n --threads {resources.cores} -o {output.bam} {input.bam}"


rule filter_bamfile:
    input:
        os.path.join(WD,'mapped_reads', 'bams','gene_counts_{sample}.name_sorted.bam'),
    output:
        os.path.join(WD,'mapped_reads','bams','{sample}_filtered.bam')
    conda:
        "envs/signature_reads.yaml"
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
        msamtools filter -b -l {params.min_len} -p {params.min_iden} -z {params.min_map} --besthit {input} > {output}
        """  

rule profile: 
    input:
        os.path.join(WD,'mapped_reads', 'bams','{sample}_filtered.bam'),
    output:
        os.path.join(WD,'signature_reads','profiles','{sample}.profile')
    conda:
        "envs/signature_reads.yaml"
    resources:
        cores = 8,
        mem_gb = 20,
        runtime = 7200 #2h in s
    params:
        multi = param_dict['multi']
    threads: 8
    shell:
        """
        msamtools profile --label={wildcards.sample} --multi={params.multi} --unit=ab \
                        --nolen -o {output}.gz {input}
        gunzip {output}.gz
    """

# Create the header for the gene count matrix
rule create_header:
    input: 
        expand(os.path.join(WD,'signature_reads','profiles','{sample}.profile'), sample=SAMPLES)
    output:
        header = temp(os.path.join(WD, 'signature_reads', 'header.txt'))
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

rule get_genes:
    input:
        os.path.join(WD,'signature_reads','profiles','{sample}.profile')
    output:
        temp(os.path.join(WD,'signature_reads','gene_names_{sample}.txt'))
    wildcard_constraints:
        sample=SAMPLES[0]
    resources:
        cores = 1,
        mem_gb = 4,
        runtime = 7200 #2h in s
    shell:
        """
        grep -v "#" {input} | awk '{{print $1}}' > {output}
        """

rule get_counts:
    input:
        os.path.join(WD,'signature_reads','profiles','{sample}.profile')
    output:
        temp(os.path.join(WD,'signature_reads','counts','{sample}.txt'))
    resources:
        cores = 1,
        mem_gb = 4,
        runtime = 7200 #2h in s
    shell:
        """
        grep -v "#" {input} | awk '{{print $2}}' > {output}
        """

rule count_matrix:
    input:
        readcounts=expand(os.path.join(WD,'signature_reads','counts','{sample}.txt'),sample=SAMPLES),
        gene_names=expand(os.path.join(WD,'signature_reads','gene_names_{sample}.txt'),sample=SAMPLES[0]),
        header = os.path.join(WD, 'signature_reads', 'header.txt')
    output:
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv')
    resources:
        cores = 1,
        mem_gb = 188,
        runtime = '3600' #1h in s
    shell:
        "paste {input.gene_names} {input.readcounts} | cat {input.header} - > {output}"
