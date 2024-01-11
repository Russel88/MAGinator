import pandas as pd
import concurrent.futures
import os

# Set the parameters from the snakefile
cov_thr = int(snakemake.params['min_reads'])
in_file = snakemake.input['cov_file']
out_file = snakemake.output['out_file']
name_file = snakemake.output['gene_names']

# Read the file into a data frame
df = pd.read_csv(in_file, sep='\t', header=0, dtype=str)

# Keep only the columns we need: "#rname", "numreads", "coverage"
df = df[["#rname", "numreads", "coverage"]]

# Calculate the count values based on the coverage threshold
df['count'] = df.apply(lambda row: int(0) if float(row['coverage']) <= cov_thr else int(row['numreads']), axis=1)

# Sort the dataframe by gene name and set the gene names as index
count_df = df.sort_values('#rname')[["#rname", "count"]]

# Save the count dataframe to a tab-separated file
count_df[["count"]].to_csv(out_file, sep='\t', index=False,header=False)
count_df[["#rname"]].to_csv(name_file, sep='\t', index=False,header=False)