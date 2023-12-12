#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2023

@author: pabloati
"""

import pandas as pd

def modify_matrix(input_file, output_file, min_reads):
    # Load the TSV file as a pandas DataFrame
    matrix = pd.read_csv(input_file, sep='\t', index_col=0)

    # Change values below the threshold to 0
    matrix[matrix < min_reads] = 0

    # Save the modified matrix to a new TSV file
    matrix.to_csv(output_file, sep='\t')

# Example usage
input_file = snakemake.input[0]
output_file = snake.make.output[0]
min_reads = snakemake.params.min_reads

modify_matrix(input_file, output_file, min_reads)

