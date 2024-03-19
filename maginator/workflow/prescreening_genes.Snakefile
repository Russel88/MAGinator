
import os

WD = config['wd']
PARAMS = config['params']
cluster_DIR = os.path.join(WD, 'signature_genes', 'clusters')

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}


with open(param_dict['reads']) as f:
    SAMPLES = ([row.split(',')[0] for row in f])

SAMPLES=set(SAMPLES)


rule all:
    input:
        os.path.join(WD, 'genes', 'representative_genes.tsv'),
        os.path.join(WD, 'signature_genes', 'clusters_sorted.RDS')
        # os.path.join(WD, 'signature_genes', 'cluster.RDS'),


# Identifying the genes that cluster across the metagenomic species / species collections
rule geneID_collectionID:
    input: 
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'genes', 'all_genes_cluster.tsv'),
        param_dict['vamb_clusters']
    output:
        os.path.join(WD, 'genes', 'representative_genes.tsv'),
        os.path.join(WD, 'genes', 'gene_lists', 'geneID_collectionID.tsv')
    conda:
        "envs/signature_genes.yaml"
    resources:
        cores = 1,
        mem_gb = 188,
        runtime = 43200 #12h in s
    script: 
        "scripts/species_collections.py"


# creating a gene count matrix only containing the genes, that does not cluster across metagenomic species / species collection
rule sort_genes_across_MGS:
    input: 
        os.path.join(WD, 'genes', 'matrix', 'gene_count_matrix.tsv'),
        os.path.join(WD, 'genes', 'representative_genes.tsv')
    output:
        os.path.join(WD, 'genes', 'matrix', 'small_gene_count_matrix.tsv')
    conda:
        "envs/signature_genes.yaml" 
    resources:
        cores = 1, 
        mem_gb = 188,
        runtime = 43200 #12h in s
    script:
        "scripts/sort_gene_mat.py"


#Converting the gene count matrix to cluster-matrices in R dataformat 
rule format_conversion:
    input:
        gene_clusters = os.path.join(WD, 'genes', 'gene_lists', 'geneID_collectionID.tsv'),
        matrix = os.path.join(WD, 'genes', 'matrix', 'small_gene_count_matrix.tsv'),
        gene_lengths = os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta.fai')
    output:
        R_clusters = os.path.join(WD, 'signature_genes', 'cluster.RDS'),
        R_gene_lengths = os.path.join(WD, 'signature_genes', 'gene_lengths.RDS'),
        clusters_dir = (directory(cluster_DIR))
    conda:
        "envs/signature_genes.yaml" 
    resources:
        cores = 1,
        mem_gb = 188,
        runtime = 86400 #1d in s
    script:
        "scripts/matrix2SG_formatconversion.R"


# sorting the genes according to the CCP
rule prescreening_genes:
    input:
        clusters = os.path.join(WD, 'signature_genes', 'cluster.RDS'),
        gene_lengths = os.path.join(WD, 'signature_genes', 'gene_lengths.RDS')
    output:
        clusters_sorted = os.path.join(WD, 'signature_genes', 'clusters_sorted.RDS'),
    params:
        min_mapped_signature_genes = param_dict['min_mapped_signature_genes'],
        n_genes = param_dict['num_signature_genes']
    conda:
        "envs/signature_genes.yaml"
    resources:
        cores = 1,
        mem_gb = 188,
        runtime = 86400 #1d in s
    script: 
        "scripts/prescreening_genes.R"
