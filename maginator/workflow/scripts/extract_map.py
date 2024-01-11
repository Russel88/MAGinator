# Import necessary modules
import os
import pysam
import sys
import multiprocessing


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
    for read in bam_file:
        # Check if the read passes the mapped percentage threshold
        if calculate_mapped_percentage(read, threshold):
            # Write the read to the filtered BAM file
            filtered_bam_file.write(read)

    # Close the input and output BAM files
    bam_file.close()
    filtered_bam_file.close()


# Extract command line arguments
bam_file_path = snakemake.input[0]
threshold = float(snakemake.params["min_map"])
output_file = snakemake.output[0]

filter_bam_file(bam_file_path, threshold, output_file)
