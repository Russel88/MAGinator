import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:   
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

rule all:
    input:
        os.path.join(WD, 'phylo', 'intermediate',  'clusterID_markerGene_geneCluster_rootGene.tsv'),
        os.path.join(WD, 'phylo', 'intermediate', 'markers_roots.fna'),
        os.path.join(WD, 'phylo', 'intermediate',  'cluster_genes.bed'),
        os.path.join(WD, 'phylo', 'intermediate',  'cluster_genes.fna.fai')

rule markers:
    input:
        os.path.join(WD, 'phylo', 'intermediate',  'gtdbtk_summary.tsv'),
        os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        os.path.join(WD, 'phylo', 'intermediate',  'gtdb_unique_bac_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate',  'gtdb_unique_ar_markers.tab'),
        os.path.join(WD, 'phylo', 'intermediate',  'gtdb_markers_bins_geneID.tsv'),
        os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta.fai')
    output:
        os.path.join(WD, 'phylo', 'intermediate',  'clusterID_markerGene_geneCluster_rootGene.tsv')
    params:
        min_gtdb_markes=param_dict['min_gtdb_markers'],
        marker_gene_cluster_prevalence=param_dict['marker_gene_cluster_prevalence']
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        mem_gb=10,
        runtime=36000 #10h in s
    script:
        "scripts/marker_genes.py"

rule fasta1:
    input:
        tab=os.path.join(WD, 'phylo', 'intermediate',  'clusterID_markerGene_geneCluster_rootGene.tsv'),
        fasta=os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta')
    output:
        os.path.join(WD, 'phylo', 'intermediate', 'markers_clusters.fna')
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        mem_gb=10,
        runtime=36000 #10h in s
    shell:
        "perl -ne 'if(/^>(\S+)/){{$c=$i{{$1}}}}$c?print:chomp;$i{{$_}}=1 if @ARGV' <(cut -f4 {input.tab} | sort | uniq) {input.fasta} > {output}"

rule fasta2:
    input:
        tab=os.path.join(WD, 'phylo', 'intermediate',  'clusterID_markerGene_geneCluster_rootGene.tsv'),
        fasta=os.path.join(WD, 'genes', 'all_genes.fna')
    output:
        os.path.join(WD, 'phylo', 'intermediate', 'markers_roots.fna')
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        mem_gb=10,
        runtime=36000 #10h in s
    shell:
        "perl -ne 'if(/^>(\S+)/){{$c=$i{{$1}}}}$c?print:chomp;$i{{$_}}=1 if @ARGV' <(cut -f5 {input.tab} | sort | uniq) {input.fasta} > {output}"

rule fasta3:
    input:
        tab=os.path.join(WD, 'tabs',  'signature_genes_cluster.tsv'),
        fasta=os.path.join(WD, 'genes', 'all_genes_nonredundant.fasta')
    output:
        os.path.join(WD, 'genes', 'signature_genes.fasta')
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        mem_gb=10,
        runtime=36000 #10h in s
    shell:
        "perl -ne 'if(/^>(\S+)/){{$c=$i{{$1}}}}$c?print:chomp;$i{{$_}}=1 if @ARGV' <(cut -f1 {input.tab}) {input.fasta} > {output}"

rule bed:
    input:
        sig=os.path.join(WD, 'tabs',  'signature_genes_cluster.tsv'),
        tab=os.path.join(WD, 'phylo', 'intermediate',  'clusterID_markerGene_geneCluster_rootGene.tsv')
    output:
        beds=os.path.join(WD, 'phylo', 'intermediate',  'signature_genes.bed'),
        bedm=os.path.join(WD, 'phylo', 'intermediate',  'marker_genes.bed')
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        mem_gb=10,
        runtime=36000 #10h in s
    shell:
        """
        awk '{{print $1"\t0\t1000000"}}' {input.sig} > {output.beds}
        cut -f4 {input.tab} | sort | uniq | awk '{{print $1"\t0\t1000000"}}' > {output.bedm}
        """

rule uniq:
    input:
        beds=os.path.join(WD, 'phylo', 'intermediate',  'signature_genes.bed'),
        bedm=os.path.join(WD, 'phylo', 'intermediate',  'marker_genes.bed'),
        fastas=os.path.join(WD, 'genes', 'signature_genes.fasta'),
        fastac=os.path.join(WD, 'phylo', 'intermediate', 'markers_clusters.fna')
    output:
        bed=os.path.join(WD, 'phylo', 'intermediate',  'cluster_genes.bed'),
        fasta=os.path.join(WD, 'phylo', 'intermediate',  'cluster_genes.fna')
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        mem_gb=10,
        runtime=36000 #10h in s
    shell:
        """
        cat {input.beds} {input.bedm} | sort | uniq > {output.bed}
        cat {input.fastas} {input.fastac} | awk '/^>/{{f=!d[$1];d[$1]=1}}f' > {output.fasta}
        """

rule index:
    input:
        os.path.join(WD, 'phylo', 'intermediate',  'cluster_genes.fna')
    output:
        os.path.join(WD, 'phylo', 'intermediate',  'cluster_genes.fna.fai')
    conda:
        "envs/phylo.yaml"
    resources:
        cores=1,
        mem_gb=10,
        runtime=36000 #10h in s
    shell:
        """
        samtools faidx {input}
        """
