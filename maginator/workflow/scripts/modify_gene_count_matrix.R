#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# Created 2023
# 
# @author: pabloati


suppressMessages(library(tidyverse))

cov.thr <- snakemake@params$min_cov
files.list <- snakemake@input[["cov_files"]]
out.file <- snakemake@output[["gene_matrix"]]


count.list <- lapply(files.list, function(file, thr = cov.thr) {
    sample.name <- gsub(".coverage", "", basename(file))
    df <- read_tsv(file, col_names = TRUE, show_col_types = FALSE) %>% 
        # Setting to 0 the counts of genes with coverage below the threshold
        mutate(!!sample.name := ifelse(coverage <= thr, 0, numreads)) %>% 
        # Genes have to be organized by name so the merging is done correctly
        arrange(`#rname`) %>%
        column_to_rownames("#rname") %>%
        select(!!sample.name)
    return(df)
})

# Merging the different samples data frames together
count.df <- Reduce(cbind,count.list) %>% rownames_to_column("Gene")
write_tsv(count.df, out.file)