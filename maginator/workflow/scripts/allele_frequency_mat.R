library(dplyr)
gene_stats_file <- snakemake@input[["gene"]]
gene_stats <- read.table(file=gene_stats_file, header=TRUE)

stats_file <- snakemake@input[["stat"]]
stats <- read.table(file=stats_file, header=TRUE)

af_cutoff <- as.numeric(snakemake@params[["af_cutoff"]])
min_signature_genes <- as.integer(snakemake@params[["min_signature_genes"]])
min_nonN <- as.numeric(snakemake@params[["min_nonN"]])
min_marker_genes <- as.integer(snakemake@params[["min_marker_genes"]])

# Identifying the samples in each tree fulfill the default tree-criteria
#(<50 detected SG, min 2 marker genes, min 50% non-n coverage)
filtered_data <- stats %>%
  filter(Marker_N >= min_marker_genes, Signature_N>=min_signature_genes, NonN > min_nonN)

# Count the number of different samples for each cluster
sample_counts <- filtered_data %>%
  group_by(Cluster) %>%
  summarize(Count = n_distinct(Sample))

# Get unique clusters
unique_clusters <- unique(gene_stats$Cluster)
unique_samples <- unique(gene_stats$Sample)

# Create an empty matrix
empty_matrix <- matrix(NA, nrow = length(unique_samples), ncol = length(unique_clusters),
                       dimnames = list(unique_samples, unique_clusters))

# Loop through each cluster
for (cluster_id in unique_clusters) {
  if ((dim(gene_stats[gene_stats$Cluster == cluster_id,]))[1]>af_cutoff){
    # Subset data for the current cluster
    cluster_data <- gene_stats[gene_stats$Cluster == cluster_id, ]
    samples <- (filtered_data[filtered_data$Cluster==cluster_id, "Sample"])
    absent_samples <- unique_samples[!unique_samples %in% samples]
    
    # Calculate the median Signature_AF per sample for the current cluster
    median_signature_af_per_sample <- cluster_data %>%
      group_by(Sample) %>%
      summarise(median_signature_af = median(Signature_AF, na.rm = TRUE))
    
    # Convert median_signature_af_per_sample to a data frame
    median_signature_af_per_sample_df <- as.data.frame(median_signature_af_per_sample)
    # Replace values based on conditions
    median_signature_af_per_sample_df$median_signature_af[median_signature_af_per_sample_df$median_signature_af > af_cutoff] <- 2
    median_signature_af_per_sample_df$median_signature_af[median_signature_af_per_sample_df$median_signature_af == af_cutoff] <- 1
    
    # Filling in the numbers in the matrix
    empty_matrix[,cluster_id] = median_signature_af_per_sample_df$median_signature_af
    empty_matrix[absent_samples,cluster_id] = 0 # if the sample is not in the tree - set NA as AF
  }
}

AF_matrix <- empty_matrix

# Save the data frame to a CSV file
write.csv(AF_matrix, file = snakemake@output[["af_matrix"]], row.names = TRUE)
