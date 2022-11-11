import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

with open(param_dict['reads']) as f:
    SAMPLES = ([row.split(',')[0] for row in f])

CLUSTERS = set(glob_wildcards(os.path.join(WD, 'signature_genes', 'clusters', 'Cluster{cluster}.RDS')).cluster)
CLUSTERS = {x for x in CLUSTERS if x.isdigit()}

rule all:
    input:
        expand(os.path.join(WD, 'signature_genes','counts', 'cluster_{cluster}_counts.RDS'), cluster=CLUSTERS),
        os.path.join(WD, 'abundance', 'abundance_phyloseq.RData')


# Identifying the refined sets of signature genes using the gene count matrix and gene lengths
rule SG_refinement:
    input:
        clusters_dir = os.path.join(WD, 'signature_genes', 'clusters', 'Cluster{cluster}.RDS'),
        clusters_sorted = os.path.join(WD, 'signature_genes', 'clusters_sorted.RDS'),
        gene_lengths = os.path.join(WD, 'signature_genes', 'gene_lengths.RDS')
    output:
        cluster_screened = os.path.join(WD, 'signature_genes', 'screened', 'cluster_{cluster}_screened.RDS'),
    conda:
        "envs/signature_genes.yaml"
    resources:
        cores = 1,
        memory = 188,
        runtime = '12:00:00'
    script: 
        "scripts/species_collections.py"


# creating a gene count matrix only containing the genes, that does not cluster across metagenomic species / species collection
rule sort_genes_across_MGS:
    input: 
        os.path.join(WD, 'genes', 'gene_count_matrix.tsv'),
        os.path.join(WD, 'genes', 'representative_genes.tsv')
    output:
        os.path.join(WD, 'genes', 'small_gene_count_matrix.tsv')
    conda:
       	"envs/signature_genes.yaml" 
    resources:
        cores = 1,
        memory = 188,
        runtime = '10:00:00'
    params:
        functions = "Functions_v4.R"
    script:
        "scripts/SG_refinement.R"


# insert rule for creating the MGS_object.RDS
rule gene_counts:
    input:
        cluster_screened = os.path.join(WD, 'signature_genes', 'screened', 'cluster_{cluster}_screened.RDS'),
        gene_lengths = os.path.join(WD, 'signature_genes', 'gene_lengths.RDS'),
        clusters_sorted = os.path.join(WD, 'signature_genes', 'clusters_sorted.RDS')
    output:
        cluster_counts = temp(os.path.join(WD, 'signature_genes','counts', 'cluster_{cluster}_counts.RDS'))
    conda:
        "envs/signature_genes.yaml"
    resources:
        cores = 1,
        memory = 188,
        runtime = '12:00:00'
    script:
        "scripts/MGS_counts.R"


# Identifying the refined sets of signature genes using the gene count matrix and gene lengths
rule SG_refinement:
    input:
        R_clusters = os.path.join(WD, 'signature_genes', 'cluster.RDS'),
        R_gene_lengths = os.path.join(WD, 'signature_genes', 'gene_lengths.RDS')
    output:
        screened_clusters = os.path.join(WD, 'signature_genes', 'clusters_screened.RDS'),
        MGS_object = os.path.join(WD, 'signature_genes', 'MGS_object.Rdata')
    conda:
       	"envs/signature_genes.yaml" 
    resources:
        cores = 8,
        memory = 188,
        runtime = '5:00:00:00' 
    params:
        functions = "Functions_v4.R"
    script:
        "scripts/SG_refinement.R"

# Creating abundance profiles from the SG
rule abundance_profile:
    input:
        R_gene_lengths = os.path.join(WD, 'signature_genes', 'gene_lengths.RDS'),
        screened_clusters = expand(os.path.join(WD, 'signature_genes', 'screened', 'cluster_{cluster}_screened.RDS'), cluster=CLUSTERS),
        R_clusters = os.path.join(WD, 'signature_genes', 'cluster.RDS'),
        annotation = os.path.join(WD, 'tabs', 'metagenomicspecies.tab')
    output:
        physeq_abundance = os.path.join(WD, 'abundance', 'abundance_phyloseq.RData'),
        tax_matrix = os.path.join(WD, 'tabs', 'tax_matrix.tsv')
    conda:
        "envs/signature_genes.yaml"
    resources:
        cores = 6,
        memory = 188,
        runtime = '12:00:00'
    script:
        "scripts/abundance_profiles.R"
