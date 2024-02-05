import pandas as pd
import os
import random
import time

# Set the parameters from the snakefile
cov_thr = int(snakemake.params['min_reads'])
in_file = snakemake.input['cov_file']
out_file = snakemake.output['out_file']
name_file = snakemake.output['names_file']
benchmark = snakemake.params['benchmark'] == "True"
map_thr = snakemake.params['min_map']
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
    count_file = "cov_pass_genes.tsv"
    count = 0
    while True:
        try:
            # Open the count file as a pandas data frame
            res_df = pd.read_csv(count_file, sep='\t',header=0,index_col=0)
            break
        except FileNotFoundError:
            # Create an empty pandas data frame, with an index named "Sample"
            res_df = pd.DataFrame()
            break
        except pd.errors.EmptyDataError:
            time.sleep(random.randint(1, 10) / 10)
            print("Waiting for the file to be released")
            count=count+1
            if count == 5:
                print("The file is locked for more than 5 iterations, please check the pipeline")
                break
    # Save the modified data frame back to the count file
    if map_filter == "pablo":
        col_name = map_thr + "map_" + str(cov_thr) + "cov"
    else:
        col_name = "Shiraz_" + str(cov_thr) + "cov"
    res_df.loc[sample_name, col_name] = num_genes
    res_df.to_csv(count_file, sep='\t', index=True)

    # while True:
    #         try:
    #             with open(count_file, 'a+') as f:
    #                 # Open the count file as a pandas data frame
    #                 fcntl.flock(f, fcntl.LOCK_EX)
    #                 res_df = pd.read_csv(count_file, sep='\t',header=0,index_col=0)
    #                 print("Writing in the file")
    #                 if filter == "pablo":
    #                     col_name = str(map_thr) + "map_" + str(cov_thr) + "cov"
    #                 else:
    #                     col_name = "Shiraz_" + str(cov_thr) + "cov"
    #                 res_df.loc[sample_name, col_name] = num_genes
    #                 res_df.to_csv(count_file, sep='\t', index=True)
    #                 fcntl.flock(f, fcntl.LOCK_UN) # release the lock
    #                 f.close()
    #                 breako
    #         except IOError as e:
    #             if e.errno == 11: # Resource temporarily unavailable (file is locked)
    #                 time.sleep(random.randint(1, 10) / 10)
    #                 print("Waiting for the file to be released")
    #         except pd.errors.EmptyDataError: # In case it is the first running, so the file can be created
    #             with open(count_file, 'a+') as f:
    #                 print("Emtpy file")
    #                 fcntl.flock(f, fcntl.LOCK_EX)
    #                 res_df = pd.DataFrame()
    #                 if filter == "pablo":
    #                     col_name = str(map_thr) + "map_" + str(cov_thr) + "cov"
    #                 else:
    #                     col_name = "Shiraz_" + str(cov_thr) + "cov"
    #                 res_df.loc[sample_name, col_name] = num_genes
    #                 res_df.to_csv(count_file, sep='\t', index=True)
    #                 fcntl.flock(f, fcntl.LOCK_UN) # release the lock
    #                 f.close()
    #                 break
        
