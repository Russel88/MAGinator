#!/usr/bin/env python3

import sys
import re

import pandas as pd

from collections import Counter
from itertools import chain

'''
Worflow of this script, after reading and preparing input data:
1. Finding possible roots for each cluster, from the taxonomy information (Optimal: same family, different genus))
2. Find common GTDB-tk marker genes for each cluster 
3.1. Filter the possible roots to only those sharing marker genes with the cluster
3.2. Pick root for each cluster (Most shared markers)
4.1 Find the geneIDs for the roots (redundant gene catalog)
4.2 Find the geneIDs for the marker genes for clusters (non-redundant gene catalog)
'''

#######################################################################################################################
# Constants
marker_gene_cluster_prevalence = float(snakemake.params.marker_gene_cluster_prevalence)
min_gtdb_markers = int(snakemake.params.min_gtdb_markes)
output_file = snakemake.output[0]

#######################################################################################################################
# Exceptions
class MarkerGeneError(Exception):
    pass

#######################################################################################################################
# Read data
tax = pd.read_csv(snakemake.input[0], sep='\t', header=None, names=('user_genome', 'classification'))

tax_cl = pd.read_csv(snakemake.input[1], sep='\t', header=None, names=('cluster', 'classification', 'clusters'))

try:
    unq_bac = pd.read_csv(snakemake.input[2], sep='\t', header=None)
    unq_bac = unq_bac[unq_bac[0].isin([y for x,y in zip(tax['classification'],tax['user_genome']) if 'd__Bacteria' in x])]
except Exception:
    unq_bac = pd.DataFrame()

try:
    unq_ar = pd.read_csv(snakemake.input[3], sep='\t', header=None)
    unq_ar = unq_ar[unq_ar[0].isin([y for x,y in zip(tax['classification'],tax['user_genome']) if 'd__Archaea' in x])]
except Exception:
    unq_ar = pd.DataFrame()

gtdb_markers_bin = pd.read_csv(snakemake.input[4], sep=r'\s+', header=None, names=('Gene','Marker','Bin','GeneCluster'), index_col=False)

gene_clusters_nr = pd.read_csv(snakemake.input[5], sep='\t', header=None, names=('Name','Length', 'Offset', 'Linebases', 'Linewidth'))

#######################################################################################################################
# Prepare data

# Make non-redundant gtdb_markers_bin file
gtdb_markers_bin_nr = gtdb_markers_bin[gtdb_markers_bin['GeneCluster'].isin(gene_clusters_nr['Name'])]

# Combine markers
unq_marker_df = pd.concat([unq_bac, unq_ar], axis=0)

# Get unique marker genes
unq_marker_dict = unq_marker_df.set_index(0)[1].to_dict()
unq_marker_dict = {k: set(v.split(',')) for (k, v) in unq_marker_dict.items()}

# Bin cluster matching
bin_clust_dict_inv = {k: v for k, v in zip(unq_marker_dict.keys(), [re.search(r'_([0-9]*)$', x).group(1) for x in unq_marker_dict.keys()])}

# Inverse above dict for better lookup
bin_clust_dict = {}
for k, v in bin_clust_dict_inv.items():
    bin_clust_dict[v] = bin_clust_dict.get(v, []) + [k]

# Split taxonomy to columns
tax_split = pd.DataFrame([x.split(';') for x in tax['classification']])

tax = pd.concat([tax.reset_index(drop=True), tax_split.reset_index(drop=True)], axis=1)

tax_clust = pd.DataFrame([x.split(';') for x in tax_cl['classification']])
tax_clust['cluster'] = tax_cl['cluster']

#######################################################################################################################
# For each cluster find possible roots
# also save domain data for later
'''
root_dict is a dict with the cluster number as key
and a nested list as value
with the first level in the list being the different taxonomic levels, family, order, class, phylum, domain
and the second level containing names of bins, which are possible roots

E.g. Assuming entire taxonomy is known:
Bins in root_dict[clust][0] have similar family as cluster, but different genus
Bins in root_dict[clust][1] have similar order as cluster, but different family
'''
root_dict = {}
domain_dict = {}
for row in tax_clust.iterrows():
    # Start level is family
    level = 4
    # Make sure taxonomy is known
    # Work up taxonomy if not known
    unknown_tax = True
    
    while unknown_tax:
        if row[1][level] is None:
            level -= 1
            if level < 0:
                sys.exit('Taxonomy unknown')
        else:
            unknown_tax = False

    # Find clusters similar at the level, but different at the more specific level (genus first iter)
    root_list = list()
    while level > -1:
        this = list(tax[(tax[level] == row[1][level]) & (tax[level+1] != row[1][level+1])]['user_genome'])
        root_list.append(this)
        level -= 1
    
    root_dict[row[1]['cluster']] = root_list
    domain_dict[row[1]['cluster']] = row[1][0]


#######################################################################################################################
# Marker genes per cluster
'''
marker_genes_dict is a dict with the cluster number as key
and a set as value
which contains the marker genes present in a fraction higher than "marker_gene_cluster_prevalence" in the cluster
'''
marker_genes_dict = {}
for clust in set(tax_cl['cluster']):

    # List of bins in cluster
    bin_list = bin_clust_dict[str(clust)] 
   
    # Get markers for each bin in cluster
    marker_list = [unq_marker_dict[x] for x in bin_list]
  
    # Only non-redundant markers for the cluster
    bin_marker_nr = gtdb_markers_bin_nr[gtdb_markers_bin_nr['Bin'].isin(bin_list)]
    marker_list = [{x for x in y if x in list(bin_marker_nr['Marker'])} for y in marker_list] 

    # Unpack and count
    marker_list_count = Counter(list(chain.from_iterable(marker_list)))

    # Get those above prevalence threshold
    marker_list_prev = [k for k, v in marker_list_count.items() if v/len(bin_list) >= marker_gene_cluster_prevalence]

    if len(marker_list_prev) == 0:
        sys.exit('No shared marker genes above prevalence threshold for cluster '+str(clust))

    marker_genes_dict[clust] = set(marker_list_prev)

#######################################################################################################################
# Find roots that are sharing marker genes
'''
final_roots is a dict with the cluster number as key
and a tuple as value
with first part containing name of bin, which should be the root (most overlapping marker genes, in most specific tax level (except genus))
and second part containing the names of the marker genes shared by the cluster and the root
'''
final_roots = dict()
for clust in set(tax_cl['cluster']):
    
    # Marker genes of focal cluster
    marker_gene_focal = marker_genes_dict[clust]

    # Marker genes of possible roots
    marker_gene_roots = [[unq_marker_dict[x] for x in level] for level in root_dict[clust]]
    bin_name_roots = [[x for x in level] for level in root_dict[clust]]

    # Get only overlapping marker genes
    marker_gene_overlap = [[list(x.intersection(marker_gene_focal)) for x in level] for level in marker_gene_roots]
  
    # At least X marker genes or max possible from focal
    min_markers = min(min_gtdb_markers, len(marker_gene_focal))

    # Only roots with at least "min_markers" overlapping genes
    marker_gene_overlap = [[(y, x) for x,y in zip(level1,level2) if len(x)>=min_markers] for level1,level2 in zip(marker_gene_overlap,bin_name_roots)]

    # Traverse up taxonomy to find the most specific root
    for level in marker_gene_overlap:
        # If no overlap check next level
        if len(level) == 0:
            continue
        # If overlap, pick the most overlapping
        else:
            ll_tmp = [len(x[1]) for x in level]
            root_bin = level[ll_tmp.index(max(ll_tmp))]        
            break
    
    final_roots[clust] = root_bin

#######################################################################################################################
# Find the geneIDs for the marker genes for the clusters and roots (bins, redundant gene catalog)
df_list = []
for clust in final_roots:

    root = final_roots[clust][0]
    markers = final_roots[clust][1]
   
    # Subset
    clust_markers = gtdb_markers_bin[(gtdb_markers_bin['Bin'].isin(bin_clust_dict[str(clust)])) & (gtdb_markers_bin['Marker'].isin(markers))]     
    bin_markers = gtdb_markers_bin[(gtdb_markers_bin['Bin'] == root) & (gtdb_markers_bin['Marker'].isin(markers))]     

    # Get most common Gene cluster and Marker gene combinations
    count_markers = Counter(list(zip(clust_markers['Marker'], clust_markers['GeneCluster'])))
    shared_markers = count_markers.most_common()
    shared_markers_df = pd.DataFrame([x[0] for x in shared_markers], columns = ['Marker','GeneCluster'])

    # Try removing all duplicates and check if enough markers are present
    try:
        shared_markers_df_copy = shared_markers_df.copy()
        shared_markers_df_copy.drop_duplicates(subset=['GeneCluster'], keep=False, inplace=True)
        shared_markers_df_copy.drop_duplicates(subset=['Marker'], keep=False, inplace=True)
        
        if len(shared_markers_df_copy) < min_gtdb_markers:
            raise MarkerGeneError('Fewer than '+str(min_gtdb_markers)+' marker genes uniquely matching gene clusters for cluster '+str(clust))
    
    # If not possible, keep first if there are duplicates
    # This means some marker genes can be represented by multiple gene clusters
    # and that gene clusters can match to multiple marker genes
    # In case of duplicates the first, most common geneCluster-markerGene combination, is chosen
    except MarkerGeneError:
        shared_markers_df_copy = shared_markers_df.copy()
        shared_markers_df_copy.drop_duplicates(subset=['GeneCluster'], inplace=True)
        shared_markers_df_copy.drop_duplicates(subset=['Marker'], inplace=True)
            
    # Subset 
    min_markers = min(min_gtdb_markers, len(shared_markers_df_copy))
    shared_markers_sub = shared_markers_df_copy.copy().iloc[:min_markers, ]
    
    # If a root bin has two of the same marker, keep one (pseudo-random) 
    bin_markers = bin_markers.drop_duplicates(subset=['Marker']).copy()
    
    # Ensure that cluster and bin have the same markers
    bin_markers = bin_markers[bin_markers['Marker'].isin(shared_markers_sub['Marker'])]
    shared_markers_sub = shared_markers_sub[shared_markers_sub['Marker'].isin(bin_markers['Marker'])]

    # Order
    shared_markers_sub.sort_values(by='Marker', inplace=True)
    bin_markers.sort_values(by='Marker', inplace=True)
    
    # Just a sanity check
    if list(shared_markers_sub['Marker']) != list(bin_markers['Marker']):
        print('Cluster: '+str(clust))
        print(shared_markers_sub)
        print(bin_markers)
        sys.exit('Oops, cluster markers and root markers do not match up')
   
    out_df = pd.DataFrame({'Cluster': len(bin_markers)*[clust],
                           'Root': len(bin_markers)*[root],
                           'Marker': list(bin_markers['Marker']),
                           'GeneCluster': list(shared_markers_sub['GeneCluster']),
                           'GeneRoot': list(bin_markers['Gene'])})
    
    df_list.append(out_df)

out_df_concat = pd.concat(df_list)
out_df_concat.to_csv(output_file, sep='\t', index=False, header=False)

