import os

WD = config['wd']
PARAMS = config['params']

with open(PARAMS, 'r') as fh:
    fl = [x.strip().split() for x in fh.readlines()]
param_dict = {x[0]: x[1] for x in fl}

VAMB = config['vamb']

rule all:
    input:
        os.path.join(WD, 'tabs', 'gene_cluster_tax_scope.tab'),
        os.path.join(WD, 'tabs', 'gene_cluster_bins.tab'),
        os.path.join(WD, 'tabs', 'synteny_clusters.tab')

rule gene_tax:
    input:
        tax=os.path.join(WD, 'tabs', 'tax_matrix.tsv'),
        mgs=os.path.join(WD, 'tabs', 'metagenomicspecies.tab'),
        cluster=os.path.join(WD, 'genes', 'all_genes_cluster.tsv'),
        vamb=VAMB
    output:
        os.path.join(WD, 'tabs', 'gene_cluster_tax_scope.tab'),
        os.path.join(WD, 'tabs', 'gene_cluster_bins.tab')
    params:
        cutoff = param_dict['tax_scope_threshold']
    conda:
        "envs/phylo.yaml"
    resources:
        cores = 1,
        mem_gb = 80,
        runtime = 7200 #2h in s
    script:
        "scripts/gene_cluster2tax.py"

rule synteny_graph:
    input:
        cluster=os.path.join(WD, 'genes', 'all_genes_cluster.tsv'),
        faa=os.path.join(WD, 'genes', 'all_genes.faa')
    output:
        graph=os.path.join(WD, 'genes', 'synteny', 'graph.tab'),
        index=os.path.join(WD, 'genes', 'synteny', 'graph_index.tab')
    conda:
        "envs/phylo.yaml"
    resources:
        cores = 1,
        mem_gb = 40,
        runtime = 36000 #10h in s
    script:
        "scripts/synteny.py"

rule synteny_mcl:
    input:
        graph=os.path.join(WD, 'genes', 'synteny', 'graph.tab'),
        index=os.path.join(WD, 'genes', 'synteny', 'graph_index.tab')
    output:
        abc=os.path.join(WD, 'genes', 'synteny', 'graph.abc'),
        mcl=os.path.join(WD, 'genes', 'synteny', 'mcl_clusters'),
        tab=os.path.join(WD, 'tabs', 'synteny_clusters.tab')
    conda:
        "envs/phylo.yaml"
    params:
        cutoff=param_dict['synteny_adj_cutoff'],
        i=param_dict['synteny_mcl_inflation']
    resources:
        cores = 40,
        mem_gb = 180,
        runtime = 86400 #1d in s
    shell:
        """
        join -1 2 -2 1 <(join -j1 <(sed 's/.*(//;s/)//;s/,//;s/\.0$//' {input.graph} | sort -k1,1) <(sort -k1,1 {input.index}) | sort -k2,2) <(sort -k1,1 {input.index}) \
                        | awk '{{print $4,$5,$3}}' > {output.abc}
        mcl <(awk '$3 > {params.cutoff}' {output.abc}) --abc -I {params.i} -te {resources.cores} -o {output.mcl}
        awk '{{gsub("$",";"NR)}}1' {output.mcl} | awk '{{gsub(/\\t/,";"NR"\\n")}}1' | sed 's/;/\t/' > {output.tab}
        """
