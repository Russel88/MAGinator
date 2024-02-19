import pandas as pd

"""
Created on 19/02/2024

@author: pabloati98
"""

# Set the parameters from the snakefile
in_file = snakemake.input['cov_file']
out_file = snakemake.output['out_file']
name_file = snakemake.output['names_file']

# Read the file into a data frame
df = pd.read_csv(in_file, sep='\t', header=0, dtype=str)

# Sort the dataframe by gene name and set the gene names as index
count_df = df.sort_values('#rname')[["#rname", "numreads"]]

# Save the count dataframe to a tab-separated file
count_df[["numreads"]].to_csv(out_file, sep='\t', index=False,header=False)
count_df[["#rname"]].to_csv(name_file, sep='\t', index=False,header=False)