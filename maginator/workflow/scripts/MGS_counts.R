Clusterlist <- readRDS(snakemake@input[["clusters_sorted"]])
GeneLengths <- readRDS(snakemake@input[["gene_lengths"]])
screened_cluster <- readRDS(snakemake@input[["cluster_screened"]])
n_genes <- as.integer(snakemake@params[["n_genes"]])

MGS_object <- list()
MGS_object[['i']] <- list()

if (length(screened_cluster) == 0){
saveRDS(MGS_object, file=snakemake@output[["cluster_counts"]])
} else{

all_genes <- which(names(GeneLengths) %in% rownames(Clusterlist[[screened_cluster$id]])) # getting all gene index numbers of the Cluster
  
for(entry in names(screened_cluster$genes)){
  sig_genes <- which(names(GeneLengths) %in% screened_cluster$"genes"[[entry]])
  MGS_object[['i']][[paste(screened_cluster$id,entry,sep='_')]] <- c(sig_genes, all_genes[!all_genes%in%sig_genes]) 
 
 if(length(all_genes)>n_genes){
  random_sig_genes <- sample(all_genes, n_genes)
} else{random_sig_genes <- sample(all_genes, length(all_genes))}
  MGS_object[['i']][[paste(screened_cluster$id,'random',sep='_')]] <-c(random_sig_genes, all_genes[!all_genes%in%random_sig_genes])
}

saveRDS(MGS_object, file=snakemake@output[["cluster_counts"]])
}
