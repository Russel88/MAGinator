set.seed(1337)
window_size <- 699
turkey_factor <- 1.2

# importing the data
Clusterlist <- readRDS(snakemake@input[["clusters"]])
GeneLengths <- readRDS(snakemake@input[["gene_lengths"]])

#Number of signature genes
n.genes <- as.integer(snakemake@params[["n_genes"]])

#minimum number of mapped genes required
n.mapped.minimum <- as.integer(snakemake@params[["min_mapped_signature_genes"]])

ids <- names(Clusterlist) #the ids of the MGS


# identifying the names of genes in all MGS the dataset
Cluster_gene_names <- c()
for(Cluster in names(Clusterlist)){
  Cluster_gene_names <- c(Cluster_gene_names, rownames(Clusterlist[[Cluster]]))
}
present_genes <- GeneLengths[names(GeneLengths) %in% Cluster_gene_names]

#preselecting genes
for(id in ids){
  
  genes_r <- Clusterlist[[id]]
  final.reads <- round(genes_r / (present_genes[rownames(genes_r)] * 10^-3))
  
  #finds the frequency of each gene. We want a ´p_sig_gene~1/n_signature_gene
  frequency_of_genes <- rowSums(final.reads>0, na.rm = T)
  
  #So, we want to choose a set of genes we can purify on that 1) are usually found in high frequency and 2) are relatively consisent with each other.
  #Select high-frequency genes, none of which are outliers by turkey's method, but a bit more sensitive
  vars <- c()
  good_set_found <- F
  if(length(frequency_of_genes)>window_size+1){
    frequency_of_genes <- sort(frequency_of_genes,decreasing = T)
    for(i in seq(1,length(frequency_of_genes)-window_size)){
      in_window <- frequency_of_genes[seq(i,i+window_size)]
      iqr <- IQR(in_window,na.rm = T)
      upper_bound <- median(in_window,na.rm = T)+iqr*turkey_factor
      lower_bound <- median(in_window,na.rm = T)-iqr*turkey_factor
      
      vars <- c(vars,var(in_window),na.rm=T)
      #fewer than 1% of samples are outliers
      if(sum(in_window>upper_bound | in_window<lower_bound)<=window_size/100){
        good_set_found <- T
        break
      }
    }
    
    if(good_set_found==F){
      names(vars) <- seq(1,length(frequency_of_genes)-window_size)
      best_vars <- na.omit(names(vars)[vars==min(vars)])
      if(length(best_vars)>0){
        i <- best_vars[1]
      }
    }
  }
  
  
  #if you can't find a good region, go on then
  if(good_set_found==T){
    acceptable_genes <- names(frequency_of_genes)[seq(i,window_size+i)]
    
    original_genes <- rownames(Clusterlist[[id]])[1:100]
    
    #adds the original ones. thy'r good for SOMETHING and reorders themnn to be the utuseq order just ot be saf
    acceptable_genes <- unique(c(original_genes,acceptable_genes))
    acceptable_genes <- rownames(Clusterlist[[id]])[rownames(Clusterlist[[id]])%in%acceptable_genes]
    #trims the set of possible genes down to acceptable ones
    Clusterlist[[id]] <- Clusterlist[[id]][acceptable_genes,]
    
  }
}

# Save the prescreened genes per cluster to a file
saveRDS(Clusterlist, file = snakemake@output[["clusters_sorted"]])
