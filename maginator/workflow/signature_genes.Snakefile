
import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}


with open(param_dict['reads']) as f:
    SAMPLES = ([row.split(',')[0] for row in f])

SAMPLES=set(SAMPLES)


rule all:
    input:
        os.path.join(WD, 'abundance', 'abundance_phyloseq.RData'),
        os.path.join(WD, 'genes', 'representative_genes.tsv')


# Identifying the genes that cluster across the metagenomic species / species collections
rule geneID_collectionID:
    input: 
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'genes', 'all_genes95_cluster.tsv'),
        param_dict['vamb_clusters']
    output:
        os.path.join(WD, 'genes', 'representative_genes.tsv'),
        os.path.join(WD, 'clusters', 'gene_lists', 'geneID_collectionID.tsv')
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
        runtime = '12:00:00'
    script:
        "scripts/sort_gene_mat.py"


#Converting the gene count matrix to cluster-matrices in R dataformat 
rule format_conversion:
    input:
        gene_clusters = os.path.join(WD, 'clusters', 'gene_lists', 'geneID_collectionID.tsv'),
        matrix = os.path.join(WD, 'genes', 'small_gene_count_matrix.tsv'),
        gene_lengths = os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta.fai')
    output:
        R_clusters = os.path.join(WD, 'signature_genes', 'cluster.RDS'),
        R_gene_lengths = os.path.join(WD, 'signature_genes', 'gene_lengths.RDS')
    conda:
       	"envs/signature_genes.yaml" 
    resources:
        cores = 1,
        memory = 188,
        runtime = '24:00:00'
    script:
        "scripts/matrix2SG_formatconversion.R"


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
        screened_clusters = os.path.join(WD, 'signature_genes', 'clusters_screened.RDS'),
        R_clusters = os.path.join(WD, 'signature_genes', 'cluster.RDS'),
        MGS_object = os.path.join(WD, 'signature_genes', 'MGS_object.Rdata'),
        annotation = os.path.join(WD, 'tabs', 'metagenomicspecies.tab')
    output:
        physeq_abundance = os.path.join(WD, 'abundance', 'abundance_phyloseq.RData'),
        init_physeq_abundance = os.path.join(WD, 'abundance', 'initial_abundance_phyloseq.RData'),
        tax_matrix = os.path.join(WD, 'tabs', 'tax_matrix.tsv')
    conda:
       	"envs/signature_genes.yaml" 
    resources:
        cores = 6,
        memory = 188,
        runtime = '12:00:00'
    script:
        "scripts/abundance_profiles.R"
