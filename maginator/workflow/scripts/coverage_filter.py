import pandas as pd
import concurrent.futures
import os

# Set the parameters from the snakefile
cov_thr = int(snakemake.params['min_reads'])
in_file = snakemake.input['cov_file']
out_file = snakemake.output['out_file']
name_file = snakemake.output['names_file']
benchmark = snakemake.params['benchmark'] == "True"
map_thr = snakemake.params['min_map'] + "map"
map_filter = snakemake.params['map_filter']

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


if benchmark:
    # Count the number of genes that are over the threshold
    num_genes = df[df['coverage'].astype(float) > int(cov_thr)]['coverage'].count()
    sample_name = os.path.basename(in_file).split(".")[0]
    count_file = "trial/maginator_readcov/cov_pass_genes.tsv"

    try: 
        # Open a predetermined file for writing the output
        res_df = pd.read_csv(count_file, sep='\t',header=0,index_col=0)
    except FileNotFoundError:
        # Create an empty pandas data frame, with an index named "Sample"
        res_df = pd.DataFrame()
    # Save the value of 'lost' in the column corresponding to the threshold and the row corresponding to the sample name
    if map_filter == "pablo":
        col_name = map_thr + "_" + str(cov_thr) + "cov"
    else:
        col_name = "Shiraz_" + str(cov_thr) + "cov"
    res_df.loc[sample_name, col_name] = num_genes
    # Save the modified data frame back to the count file
    res_df.to_csv(count_file, sep='\t', index=True)