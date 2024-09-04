# Loading relevant data
library(ggplot2)

GeneLengths <- readRDS(snakemake@input[["gene_lengths"]])
Clusterlist <- readRDS(snakemake@input[["clusters_sorted"]])
sg_files <- snakemake@input[["screened_clusters"]]
screened_clusters <- do.call("rbind", lapply(sg_files, readRDS))
taxmat <- read.csv(snakemake@input[["tax_matrix"]], sep="\t", header=FALSE) 

# Initializing relevant parameters
n_signature_genes_expected <- '95'
minimum_sampels <- as.integer(snakemake@params[["min_samples"]])
n.genes <- as.integer(snakemake@params[["n_genes"]])

print(snakemake@params)

#loading colors
cmcol1 = c('#0571B0', '#92C5DE', '#999999', '#E0E0E0', '#CA0020',
           '#FF8784', '#FFD900', '#FFE967', '#349F2C', '#BAE4B3',
           '#6B3D9A', '#CAB2D6')

# The function for visualizing the reads and the expected distribution along with it's MSE and taxonomy
pois.two.plot <- function(i.Genes, i.Reads, i.mse, r.Genes, r.Reads, r.mse, id){ 
  # Plotting the reads and mapped genes for each sample
  #
  # Args:
  #   Genes: A named vector containing the number of genes mapped in the samples
  #      i.Genes: initial gene names
  #      r.Genes: Refined gene names
  #   Reads: Vector containing the number of reads mapping to the Genes. Genes and Reads must 
  #           have the same length, greater than one, with no missing values
  #   mse:  The MSE of the MGS when compared to the expected distribution
  #   id is the ID of the species in the HGMGS object
  #
  # Returns:
  #   A plot of all samples with their corresponding genes and mapped reads (log)
  
  
  df_plot <- do.call(rbind, Map(data.frame, Reads = i.Reads, Genes = i.Genes, name = names(i.Genes)))
  df_r_plot <- do.call(rbind, Map(data.frame, r.Reads = r.Reads, r.Genes = r.Genes, name = names(r.Genes)))
  
  
  p = ggplot(log = "x", xlim = c(0.1, 1e5), cex = 0.8, pch = 25, ylim = c(0, n.genes), 
             data = df_plot, aes(text = name, y = Genes, x = Reads)) + 
    geom_point(size = 0.7, aes(color = "Initial")) +
    stat_function(fun =  function(Genes) (1 - ((n.genes - 1) / n.genes)^Genes) * n.genes, aes(Genes, color="Expected"), size=0.25, alpha=0.5) + #t he expected distribution
    xlab('Reads mapped (log)') +
    ylab('Signature genes detected') +
    ggtitle(id) +
    geom_point(data=df_r_plot, aes(x=r.Reads, y=r.Genes,color="Refined"), size=0.7)
  
  plotly.data <- (p + scale_x_log10()+ scale_colour_manual(values = c("Initial" = cmcol1[2], "Refined"=cmcol1[5], "Expected"="black")))
  plotly.data <- plotly.data + annotate("text",x=1,y=100,hjust=.2, label = paste("initial MSE:", round(i.mse, 2), "\n refined MSE:", round(r.mse, 2))) 
  plotly.data <- plotly.data  + theme_minimal() + theme(legend.position="bottom")
  return(plotly.data)  
}


Cluster_gene_names <- c()
for(Cluster in names(Clusterlist)){
  Cluster_gene_names <- c(Cluster_gene_names, rownames(Clusterlist[[Cluster]]))
}

present_genes <- GeneLengths[names(GeneLengths) %in% Cluster_gene_names]


pdf(snakemake@output[["plot_pdf"]])

for(Cluster in taxmat$V1){

if (length(screened_clusters[screened_clusters[,'id']==Cluster,]) > 1){

tax=taxmat[taxmat$V1==Cluster,8]
if (is.na(taxmat[taxmat$V1==Cluster,8])){
tax=taxmat[taxmat$V1==Cluster,7]}

colsum <- colSums(Clusterlist[[Cluster]][1:n.genes, ])
mapped <- colsum[colsum >= minimum_sampels]


refined_reads <- round(Clusterlist[[Cluster]][screened_clusters[screened_clusters[,'id']==Cluster,]$genes$best, names(mapped)] / 
                 +                   (present_genes[rownames(Clusterlist[[Cluster]][screened_clusters[screened_clusters[,'id']==Cluster,]$genes$best, ])] * 10^-3))


refined_countreads <- colSums(refined_reads)
refined_genes <- colSums(refined_reads > 0) 

reads <- round(Clusterlist[[Cluster]][screened_clusters[screened_clusters[,'id']==Cluster,]$genes$init,names(mapped)] / 
                 +                   (present_genes[rownames(Clusterlist[[Cluster]][screened_clusters[screened_clusters[,'id']==Cluster,]$genes$init, ])] * 10^-3))

countreads <- colSums(reads)
genes <- colSums(reads > 0) 

p <- pois.two.plot(genes, countreads, as.numeric(screened_clusters[screened_clusters[,'id']==Cluster,]$MSE$initial), refined_genes, refined_countreads, as.numeric(screened_clusters[screened_clusters[,'id']==Cluster,]$MSE$best), paste(Cluster, tax))

print(p)

}}
dev.off()
