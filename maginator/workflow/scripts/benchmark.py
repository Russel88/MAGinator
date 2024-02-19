import pandas as pd

# Set the parameters from the snakefile
in_file = snakemake.input[0]
reads_file = snakemake.output['read_counts']
genes_file = snakemake.output['gene_counts']
length = snakemake.params['min_len']
identity = snakemake.params['min_iden']
map = snakemake.params['min_map']

row_name = length + "len_" + identity + "iden_" + map + "map"

# Function to save the data
def save_data(file, data):
    # Open the reads file as a pandas data frame, and see if it is empty
    try:
        # Open the count file as a pandas data frame
        res_df = pd.read_csv(file, sep='\t',header=0,index_col=0)
    except FileNotFoundError:
        # Create an empty pandas data frame, with an index named "Sample"
        res_df = pd.DataFrame(columns=data.index)
        
    # Add the benchmark information to the data frame
    res_df = res_df._append(pd.Series(data,name=row_name))
    
    res_df.to_csv(file, sep='\t', index=True)
    return()

# Read the file into a data frame
df = pd.read_csv(in_file, sep='\t', header=0,index_col="Gene")

# Calculate the number of reads that are mapped per sample and the genes with reads
reads_counts = df.sum(axis=0)
genes_pass = (df != 0).sum()

# Save the data to the respective files

save_data(reads_file, reads_counts)
save_data(genes_file, genes_pass)