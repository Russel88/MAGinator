library(matrixStats) # colMedians

genemat <- read.table(snakemake@input[["matrix"]], header = T, sep = "\t", row.names=1)
rows <- length(genemat[,1]) #    
cols <- length(genemat[1,])

gene_clusters <- read.table(snakemake@input[["gene_clusters"]], header = F, sep="\t")
clustercounts <- as.data.frame(table(gene_clusters[,2])) #counting the genes of all bins
cluster_order <- unique(gene_clusters[,2])


colnames(gene_clusters) <- c("SEQUENCEID", "CLUSTERID")
cluster_names <- c()
dir.create(snakemake@output[["clusters_dir"]])

for (i in cluster_order){
  gene_names <- gene_clusters[gene_clusters$CLUSTERID == i,1] # getting the gene names of the cluster
  
  # Sort the gene count matrix of the MSP according to correlation
  AbuMatrix <- as.matrix(genemat[gene_names,, drop=FALSE])
  Median_profile <- colMedians(AbuMatrix)
  pcc2MedianProfile <-cor(t(AbuMatrix) , Median_profile )
  gene_order <- row.names(pcc2MedianProfile)[order(pcc2MedianProfile, decreasing = TRUE)]
  ordered_AbuMatrix <- AbuMatrix[gene_order,, drop=FALSE]
  
  m_name <- paste('Cluster', i, sep = "")
  cluster_mat <- paste('Cluster', i, sep = "")
  cluster_mat <- ordered_AbuMatrix
  assign(m_name, ordered_AbuMatrix)
  cluster_names <- c(cluster_names, m_name)

  outfile <- paste(snakemake@output[["clusters_dir"]], '/', m_name, '.RDS', sep="")
  saveRDS(cluster_mat, file=outfile)
}


 
#combining the clusters into a large list
clusterlist <- mget(cluster_names)

saveRDS(clusterlist, file=snakemake@output[["R_clusters"]])

# Saving the gene lengths as a named array of integers
gene_lengths <- read.table(snakemake@input[["gene_lengths"]], header=F, sep="\t")
genes <- split(gene_lengths[,2]*3, gene_lengths[,1]) #multiplying by 3 to go from prot to nucleotide
genes <- unlist(genes)

saveRDS(genes, file=snakemake@output[["R_gene_lengths"]])
