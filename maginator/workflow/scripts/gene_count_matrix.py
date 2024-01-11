import pandas as pd
import concurrent.futures
import os

# Set the parameters from the snakefile
cov_thr = int(snakemake.params['min_reads'])
files_list = snakemake.input['cov_files']
out_file = snakemake.output['gene_matrix']

# Function to process each file
def process_file(file):
    # Extract the sample name from the file path
    sample_name = file.split('/')[-1].replace('.coverage', '')
    # Read the file into a dataframe
    df = pd.read_csv(file, sep='\t', header=0, dtype=str)
    
    # Calculate the count values based on the coverage threshold
    df[sample_name] = df.apply(lambda row: int(0) if float(row['coverage']) <= cov_thr else int(row['numreads']), axis=1)
    
    # Sort the dataframe by gene name and set the gene name as the index
    df = df.sort_values('#rname').set_index('#rname')[[sample_name]]
    
    return df

# Initialize an empty list to store the count dataframes
count_list = []

# Create a ThreadPoolExecutor with the number of available cores
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Submit each file for processing
    futures = [executor.submit(process_file, file) for file in files_list]
    
    # Retrieve the results as they become available
    for future in concurrent.futures.as_completed(futures):
        count_list.append(future.result())

# Concatenate the count dataframes along the columns, reset the index, and rename the 'index' column to 'Gene'
count_df = pd.concat(count_list, axis=1).reset_index().rename(columns={'index': 'Gene'})

# Rename the index to 'Gene'
count_df.index.name = 'Gene'
# Save the count dataframe to a tab-separated file
count_df.to_csv(out_file, sep='\t', index=False)
