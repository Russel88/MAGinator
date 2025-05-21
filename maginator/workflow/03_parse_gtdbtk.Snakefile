import os
import distutils.util

WD = config['wd']
PARAMS = config['params']
VAMB = config['vamb']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers_bins_geneID.tsv'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdbtk_summary.tsv')

rule split_gtdbtk:
    input:
        os.path.join(WD, 'gtdbtk', 'gtdbtk.bac120.summary.tsv')
    output:
        os.path.join(WD, 'gtdbtk', 'done.txt')
    resources:
        cores=1,
        mem_gb=32,
        runtime=36000
    run:
        from collections import defaultdict

        input_file = input[0]
        done_file = output[0]
        outdir = os.path.dirname(done_file)

        with open(input_file, 'r') as f:
            lines = f.readlines()

        header = lines[0].rstrip('\n')
        suffix_to_lines = defaultdict(list)

        for line in lines[1:]:
            line = line.rstrip('\n')
            if not line.strip():
                continue
            genome = line.split('\t')[0]
            suffix = genome.split('_')[-1]
            suffix_to_lines[suffix].append(line)

        for suffix, grouped_lines in suffix_to_lines.items():
            dir_path = os.path.join(outdir, suffix)
            os.makedirs(dir_path, exist_ok=True)
            out_path = os.path.join(dir_path, 'gtdbtk.bac120.summary.tsv')
            with open(out_path, 'w') as out_f:
                out_f.write(header + '\n')
                for l in grouped_lines:
                    out_f.write(l + '\n')

        with open(done_file, 'w') as f:
            f.write("done\n")


rule parse_gtdbtk:
    input:
        os.path.join(WD, 'gtdbtk'),
        ancient(WD),
        ancient(os.path.join(WD, 'gtdbtk','done.txt'))
    output:
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'genes', 'all_genes.fna'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_unique_bac_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_unique_ar_markers.tab'),
        os.path.join(WD, 'genes', 'all_genes.faa')
    params:
        param_dict['annotation_prevalence'],
        VAMB,
        distutils.util.strtobool(param_dict['mgs_collections'])
    conda:
        "envs/filter_gtdbtk.yaml"
    resources:
        cores=1,
        mem_gb=32,
        runtime=36000 #10h in s
    script:
        "scripts/parse_gtdbtk.py"

# identifying the sequence type which will be the basis for the gene clustering
if param_dict["clustering_type"] == "protein":
    sequence_type_cluster = os.path.join(WD, 'genes', 'all_genes.faa')
elif param_dict["clustering_type"] == "nucleotide":
    sequence_type_cluster = os.path.join(WD, 'genes', 'all_genes.fna')


# Get representative genes from all genes.
rule repres_genes:
    input:
        sequence_type_cluster
    output:
        tsv = os.path.join(WD, 'genes', 'all_genes_cluster.tsv')
    resources:
        cores = 14,
        mem_gb = 50,
        runtime = 172800 #2d in s 
    params:
        tmp_dir = os.path.join(WD, 'tmp'),
        out_prefix = os.path.join(WD, 'genes', 'all_genes'),
        cov = param_dict["clustering_coverage"],
        seq_id = param_dict["clustering_min_seq_id"]
    conda:
        "envs/filter_gtdbtk.yaml"
    shell:
        "mmseqs easy-linclust --min-seq-id {params.seq_id} -c {params.cov} --threads {threads} {input} {params.out_prefix} {params.tmp_dir}; rm -r {params.tmp_dir};"


# Add gene clusters to GTDB-tk data
rule join:
    input:
        gtdb = os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers.tab'),
        cluster = os.path.join(WD, 'genes', 'all_genes_cluster.tsv')
    output:
        os.path.join(WD, 'phylo', 'intermediate', 'gtdb_markers_bins_geneID.tsv')
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 86400 #1d in s
    shell:
        "join -1 1 -2 2 <(sort -k1,1 {input.gtdb}) <(sort -k2,2 {input.cluster}) > {output}"

# Collect GTDB-tk summary info
rule collect:
    input:
        os.path.join(WD, 'gtdbtk','gtdbtk.bac120.summary.tsv')
    output:
        os.path.join(WD, 'phylo', 'intermediate', 'gtdbtk_summary.tsv')
    resources:
        cores = 1,
        mem_gb = 20,
        runtime = 86400 #1d in s
    shell:
        "cat {input} | cut -f1,2 > {output}"
