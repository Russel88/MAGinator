import os
import sys
import re
import shutil
import glob

import pandas as pd

from Bio import SeqIO
from collections import Counter

'''
Read GTDB-tk taxonomy and define metagenomic species
Collect all genes
'''

# Make directories
try:
    os.mkdir(os.path.join(snakemake.input[1], 'genes'))
except FileExistsError:
    pass
try:
    os.mkdir(os.path.join(snakemake.input[1], 'tabs'))
except FileExistsError:
    pass

def most_common(ll):
    return(Counter(ll).most_common(1)[0])

empty_tax = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

# Get consensus annotation for clusters
tax_list = []
gene_path_list = []
gene_path_prot_list = []
for clust in os.listdir(snakemake.input[0]):

    # Try reading both bacterial and archaeal summaries.
    try:
        tax_bac = pd.read_csv(os.path.join(snakemake.input[0], clust, 'gtdbtk.bac120.summary.tsv'), sep='\t', header=0)  
    except FileNotFoundError:
        tax_bac = None
    try: 
        tax_ar = pd.read_csv(os.path.join(snakemake.input[0], clust, 'gtdbtk.ar122.summary.tsv'), sep='\t', header=0)  
    except FileNotFoundError:
        tax_ar = None

    # Combine
    try:
        tax = pd.concat([tax_bac, tax_ar])
    except Exception:
        tax = None

    if tax is not None:
        
        classification = [x.split(';') for x in tax['classification']]
        
        # Remove unclassified
        classification = [x for x in classification if len(x) == 7]
       
        # Traverse from species annotation and up
        # Pick the annotation if the most common annotation is above prevalence set by parameter
        level = 6
        while level >= 0:
            annot = [x[level] for x in classification]
            cons = most_common(annot)
            if cons[1] >= float(snakemake.params[0])*len(annot):
                rep_annot = classification[annot.index(cons[0])]
                # If annotation not species level, add empty annotations for consistency
                while level < 6:
                    rep_annot[level+1] = empty_tax[level+1]
                    level += 1
                final_annot = ';'.join(rep_annot)
                break
            else:
                level -= 1

        tax_list.append([[clust], final_annot])

    # Collect gene sequence paths
    gene_path = os.path.join(snakemake.input[0], clust, 'identify', 'intermediate_results', 'marker_genes')

    for bin_id in os.listdir(gene_path):
        gene_path_list.append(os.path.join(gene_path, bin_id, bin_id+'_protein.fna'))
        gene_path_prot_list.append(os.path.join(gene_path, bin_id, bin_id+'_protein.faa'))

# If similar at species level, then aggregate
species_list = []
mgs_list = []
for clust in tax_list:
    # vamb_cluster parameter is set, not clusters are aggregated
    if bool(snakemake.params[2]):
        mgs_list.append(clust)
    # If vamb_cluster parameter is not set, aggregate if species annotation is similar
    else:
        # If not species annotation, just add to final list
        if bool(re.search('s__$', clust[1])):
            mgs_list.append(clust)
        # If species annotation, check if it's already added to list - if False add to final list. If True append to existing
        else:
            if clust[1] in species_list:
                mgs_list[[x[1] for x in mgs_list].index(clust[1])][0] += clust[0]
            else:
                mgs_list.append(clust)
                species_list.append(clust[1])

# Write metagenomicspecies file
with open(snakemake.output[0], 'w') as fh:
    for mgs in mgs_list:
        rep = min([int(x) for x in mgs[0]])
        fh.write('{}\t{}\t{}\n'.format(str(rep), mgs[1], ','.join(mgs[0])))

# Write genes
with open(snakemake.output[1], 'wb') as wfh:
    for f in gene_path_list:
        with open(f, 'rb') as rfh:
            shutil.copyfileobj(rfh, wfh)

#write genes as proteins
with open(snakemake.output[5], 'wb') as wfh:
    for f in gene_path_prot_list:
        with open(f, 'rb') as rfh:
            shutil.copyfileobj(rfh, wfh)


# Collect marker gene info for phylogenetic analyses
## Get contig-bin info
with open(snakemake.params[1], 'r') as fh:
    vamb = [x.strip().split() for x in fh.readlines()]
vamb = {x[1]: x[0] for x in vamb}

## Find all GTDB marker files
with open(snakemake.output[2], 'w') as wfh:
    for f in glob.glob(os.path.join(snakemake.input[0], '*', 'identify', 'intermediate_results', 'marker_genes', '*', '*tophit.tsv')):
        with open(f, 'r') as rfh:
            nl = 0
            for ll in rfh:
                if nl > 0:
                    line = ll.strip().split()
                    wfh.write(line[0]+'\t'+re.sub(',.*', '', line[1])+'\t'+vamb[re.sub('_[0-9]*$', '', line[0])]+'\n')
                nl += 1

# Collect all unique markers for phylogenetic analyses
def unique_gtdb(domain, output):
    with open(output, 'w') as wfh:
        for f in glob.glob(os.path.join(snakemake.input[0], '*', 'identify', 'gtdbtk.'+domain+'.markers_summary.tsv')):
            with open(f, 'r') as rfh:
                nl = 0
                for ll in rfh:
                    if nl > 0:
                        line = ll.strip().split()
                        wfh.write(line[0]+'\t'+line[5]+'\n')
                    nl += 1

unique_gtdb('bac120', snakemake.output[3])
unique_gtdb('ar122', snakemake.output[4])
