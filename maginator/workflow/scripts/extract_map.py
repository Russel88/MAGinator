# Import necessary modules
import os
import pysam
import sys
import multiprocessing
import pandas as pd

# Define a determine if the read is over the percentage of mapped bases
def calculate_mapped_percentage(read, threshold):
    # Calculate the percentage of mapped bases
    mapped_bases = read.query_alignment_length
    total_bases = read.query_length
    mapped_percentage = (mapped_bases / total_bases)
    
    # Check if the mapped percentage is greater than or equal to the threshold
    return mapped_percentage >= threshold/100

# Define a function to filter a BAM file based on the mapped percentage threshold
def filter_bam_file(bam_file_path, threshold, output_file):
    # Open the input BAM file for reading
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    # Create a new BAM file for writing the filtered reads
    filtered_bam_file = pysam.AlignmentFile(output_file, "wb", header=bam_file.header)

    # Iterate over each read in the input BAM file
    count = 0
    for read in bam_file:
        # Check if the read passes the mapped percentage threshold
        if calculate_mapped_percentage(read, threshold):
            # Write the read to the filtered BAM file
            filtered_bam_file.write(read)
        else:
            count += 1
    # Close the input and output BAM files
    bam_file.close()
    filtered_bam_file.close()
    return count



# Extract command line arguments
bam_file_path = snakemake.input[0]
threshold= float(snakemake.params["min_map"])
cov_thr = snakemake.params["min_cov"] + "cov"
output_file = snakemake.output[0]
benchmark = snakemake.params["benchmark"] == "True"

lost = filter_bam_file(bam_file_path, threshold, output_file)

if benchmark:
# Open a predetermined file for writing the output
    count_file = "map_deleted_reads.tsv"
    map_thr = str(threshold) + "map"
    sample_name = os.path.basename(bam_file_path).split(".")[0].split("_")[2]

    try:
        # Open the count file as a pandas data frame
        count_df = pd.read_csv(count_file, sep='\t',header=0,index_col=0)
    except FileNotFoundError:
        # Create an empty pandas data frame, with an index named "Sample"
        count_df = pd.DataFrame()
    # Save the value of 'lost' in the column corresponding to the threshold and the row corresponding to the sample name
    count_df.loc[sample_name, map_thr + "_" + cov_thr] = lost
    # Save the modified data frame back to the count file
    count_df.to_csv(count_file, sep='\t', index=True)    