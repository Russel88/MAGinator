#!/usr/bin/env python3

import re
import sys
import itertools

import pandas as pd
import numpy as np

####### Read #######
tax_df = pd.read_csv(snakemake.input['tax'], header=None, names=('MGS', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep='\t')
mgs_df = pd.read_csv(snakemake.input['mgs'], header=None, names=('MGS', 'Taxonomy', 'Clusters'), sep='\t')
clust_df = pd.read_csv(snakemake.input['cluster'], header=None, names=('GeneCluster', 'Gene'), sep='\t')
vamb_df = pd.read_csv(snakemake.input['vamb'], header=None, names=('Bin', 'Contig'), sep='\t')

cutoff = float(snakemake.params['cutoff'])

###### Transform ######
# Expand mgs data.frame to contain a row per cluster for easy merging
mgs_zip = zip(list(mgs_df['MGS']), list(mgs_df['Clusters']))
mgs_list = list(itertools.chain.from_iterable([[('Cluster'+str(x[0]), y) for y in str(x[1]).split(',')] for x in mgs_zip]))
mgs_df_exp = pd.DataFrame(mgs_list, columns=('MGS', 'Cluster'))

####### Merging ######
# Merge taxonomy
tax_df_new = tax_df.merge(mgs_df_exp, on='MGS', how='left')

# Fill cluster info for NAs
tax_df_new.loc[tax_df_new.Cluster.isna(), 'Cluster'] = [re.sub('Cluster', '', x) for x in list(tax_df_new.loc[tax_df_new.Cluster.isna(), 'MGS'])]

# Add contig column to gene clustering
clust_df['Contig'] = [re.sub('_[0-9]*$', '', x) for x in clust_df['Gene']]

# First merge
merge_df1 = clust_df.merge(vamb_df, on='Contig')

# Add cluster column to bin table
merge_df1['Cluster'] = [re.sub('.*_', '', x) for x in merge_df1['Bin']]

# Second merge
merge_df2 = merge_df1.merge(tax_df_new, on='Cluster')

# Write first table
merge_df1[['GeneCluster', 'Bin']].to_csv(snakemake.output[1], sep='\t', index=False)

####### Summarising ######
# Split
merge_single = merge_df2.drop_duplicates(subset='GeneCluster', keep=False)
merge_dup = merge_df2.loc[merge_df2.duplicated(subset='GeneCluster', keep=False), :]

count_df = pd.DataFrame(merge_dup.groupby('GeneCluster')[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']].value_counts(normalize=True,dropna=False))
count_df.reset_index(inplace=True)
count_df = count_df.rename(columns={0: 'proportion'})

# All those above cutoff are consistent at Species level
species_df = count_df[(count_df["proportion"] > cutoff) & ~(count_df['Species'].isna())]

# Traverse through taxonomy, aggregate and check if cutoff is met
# Genus level
tmp_df = count_df[~count_df['GeneCluster'].isin(species_df['GeneCluster'])].groupby(['GeneCluster', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus'], dropna=False)["proportion"].sum().reset_index()
genus_df = tmp_df[(tmp_df["proportion"] > cutoff) & ~(tmp_df['Genus'].isna())]

# Family level
tmp_df = tmp_df[~tmp_df['GeneCluster'].isin(genus_df['GeneCluster'])].groupby(['GeneCluster', 'Domain', 'Phylum', 'Class', 'Order', 'Family'], dropna=False)["proportion"].sum().reset_index()
family_df = tmp_df[(tmp_df['proportion'] > cutoff) & ~(tmp_df['Family'].isna())]

# Order level
tmp_df = tmp_df[~tmp_df['GeneCluster'].isin(family_df['GeneCluster'])].groupby(['GeneCluster', 'Domain', 'Phylum', 'Class', 'Order'], dropna=False)["proportion"].sum().reset_index()
order_df = tmp_df[(tmp_df['proportion'] > cutoff) & ~(tmp_df['Order'].isna())]

# Class level
tmp_df = tmp_df[~tmp_df['GeneCluster'].isin(order_df['GeneCluster'])].groupby(['GeneCluster', 'Domain', 'Phylum', 'Class'], dropna=False)["proportion"].sum().reset_index()
class_df = tmp_df[(tmp_df['proportion'] > cutoff) & ~(tmp_df['Class'].isna())]

# Phylum level
tmp_df = tmp_df[~tmp_df['GeneCluster'].isin(class_df['GeneCluster'])].groupby(['GeneCluster', 'Domain', 'Phylum'], dropna=False)["proportion"].sum().reset_index()
phylum_df = tmp_df[(tmp_df['proportion'] > cutoff) & ~(tmp_df['Phylum'].isna())]

# Domain level
tmp_df = tmp_df[~tmp_df['GeneCluster'].isin(phylum_df['GeneCluster'])].groupby(['GeneCluster', 'Domain'], dropna=False)["proportion"].sum().reset_index()
domain_df = tmp_df[(tmp_df['proportion'] > cutoff) & ~(tmp_df['Domain'].isna())]

# No level
tmp_df = tmp_df[~tmp_df['GeneCluster'].isin(domain_df['GeneCluster'])].groupby(['GeneCluster'])["proportion"].sum().reset_index()

summary_df = pd.concat([species_df, genus_df, family_df, order_df, class_df, phylum_df, domain_df, tmp_df])
del summary_df["proportion"]

summary_df.to_csv(snakemake.output[0], sep='\t', index=False, na_rep='NA')

