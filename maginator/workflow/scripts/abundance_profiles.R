# Converting the read count matrix to abundance profiles with GTDB-tk annotation

# Initializing 
library(phyloseq)
library(MASS)
library(stringr)

# Loading relevant input
GeneLengths <- readRDS(snakemake@input[["R_gene_lengths"]]) # The gene lengths
sg_files <- snakemake@input[["screened_clusters"]]
screened_clusters <- do.call("rbind", lapply(sg_files, readRDS))
#load(snakemake@input[["MGS_object"]]) # contain the SGs of the Clusters
Clusterlist <- readRDS(snakemake@input[["R_clusters"]]) # read count mat ofclusters
taxonomy <- read.csv(snakemake@input[["annotation"]], header=FALSE, sep="\t") # the taxonomy 
stat <- snakemake@params[["stat"]] # the statistic used to calculate the abundance
pctg <- as.integer(snakemake@params[["percentage"]])/100 # percentage of the SG distribution that will be left out to calculate abundances
colnames(taxonomy) <- c("Cluster","Taxonomy")

#setting important variables
gene_index <- seq(1,length(GeneLengths))
gene_names <- names(GeneLengths)
n.mapped.minimum <- as.integer(snakemake@params[["min_genes"]]) #The number of reads that needs reads that map to count the cluster as present
n.genes <- as.integer(snakemake@params[["n_genes"]]) # number of signature genes

# inserting NA for the Clusters that do not have a annotation
taxmat <- matrix("NA", nrow = length(names(Clusterlist)), ncol = 7)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(taxmat) <- names(Clusterlist)
MGSids <- unlist(lapply(str_split(taxonomy$Cluster,","), '[[',1))

for (cluster in names(Clusterlist)){
  cluster_no <- str_split(cluster, "Cluster")[[1]][2]

  if (cluster_no %in% MGSids){ # if the cluster is annotated
    tax <- taxonomy$Taxonomy[lapply(str_split(as.character(taxonomy$Cluster),","), '[[',1)==cluster_no]
    tax <- substr(strsplit(tax, ";")[[1]], 4, nchar(strsplit(tax, ";")[[1]]))
    row <- rep("NA", 7)
    row[1:length(tax)] <- as.character(tax)
    taxmat[cluster,] <- row
  } 
}

write.table(taxmat, file=snakemake@output[["tax_matrix"]], row.names=TRUE, col.names=FALSE, sep="\t",quote = FALSE)

# Identifying the "Maintax"
taxonomy$MainTax <- rep("NA", length(taxonomy[,1]))
count <- 0 
for (MGS in taxonomy$Cluster){
  if (length(taxonomy$Taxonomy[taxonomy$Cluster==MGS])==0){
    na_df <- data.frame(MGS, "NA;NA;NA", "NA")
    colnames(na_df) <- c("Cluster", "Taxonomy", "MainTax")
    taxonomy <- rbind(taxonomy, na_df)
    count <- count + 1
  } else {
    main_tax <- paste(tail(strsplit(as.character(taxonomy$Taxonomy[taxonomy$Cluster==MGS]),";")[[1]], n=2)[1], tail(strsplit(as.character(taxonomy$Taxonomy[taxonomy$Cluster==MGS]),";")[[1]], n=2)[2])
    taxonomy$MainTax[taxonomy$Cluster == MGS] <- main_tax
  }
}
rownames(taxonomy) <- taxonomy$Cluster # set the clusterID as rownames (for the phyloseq)
taxonomy$Cluster <- NULL #removing the column witht the clusterID


# Prepare a table containing the signature geneID's in column1 and the corresponding cluster in column2
sg_cluster <- matrix("NA", nrow = length(names(Clusterlist))*100, ncol = 2)
sg_cluster[,2]<- sort(rep(names(Clusterlist), 100))

# The read counts are normalized according to gene lengths
final.read.matrix <- matrix(NA, nrow=dim(Clusterlist[[1]])[2], ncol=length(Clusterlist))
sample.names <- colnames(Clusterlist[[1]]) 
rownames(final.read.matrix) <- sample.names
colnames(final.read.matrix) <- names(Clusterlist) 

final.Clusterlist <- Clusterlist
sg_reads <- list()
for (id in names(Clusterlist)){  
  # Repeat for the final SG
  final.gene.names <- screened_clusters[,3][screened_clusters[,'id']==id][1][[1]]$best
  final.colsum <- colSums(Clusterlist[[id]][final.gene.names, ])
  final.mapped <- names(final.colsum[final.colsum >= n.mapped.minimum])
  final.not_mapped <- sample.names[!sample.names %in% final.mapped]
  final.Clusterlist[[id]][,final.not_mapped] <- 0 # setting the counts to 0 if less than n.mapped.minimum genes have reads that map
  
  # The readcounts are divided by the gene length
  final.reads <- final.Clusterlist[[id]][final.gene.names, ] / GeneLengths[final.gene.names]
  sg_reads[[id]] <- final.reads
  # summing the read counts for the id/cluster/MGS
  if (stat == "sum"){
    abundance <- colSums(final.reads)
  } else if (stat == "tt_trun"){ # Obtain the truncated average of the read counts
    # Calculate truncated mean for each column
    abundance <- apply(final.reads, 2, function(x) {
      if (all(x == 0)) { # If all values are 0, the truncated mean would return NaN
        return(0)
      } else {
       quantiles <- quantile(x, probs = c(pctg, 1-pctg))
        return(mean(x[x >= quantiles[1] & x <= quantiles[2]])) # Should it be also equal??
      }
    })
  } else if (stat == "ot_trun"){ #One tailed truncated mean (only truncated in the top X genes)
    # Calculate truncated mean for each column
    abundance <- apply(final.reads, 2, function(x) {
      if (all(x == 0)) { # If all values are 0, the truncated mean would return NaN
        return(0)
      } else {
        quantiles <- quantile(x, probs = c(1-pctg))  # Should it be also equal??
        return(mean(x[x <= quantiles]))
      }
    })

  }
  final.read.matrix[, id] <- abundance
  if (length(final.gene.names)>0){
  if (length(final.gene.names)!=n.genes){
  final.gene.names<-c(final.gene.names, rep("NA", (100-length(final.gene.names))))}
  sg_cluster[sg_cluster[,2]==id] <- matrix(c(final.gene.names, rep(id, length(final.gene.names))), ncol=2)
}}

write.csv(final.read.matrix,"Absolute_counts.tsv",sep="\t")

final.abundance <- final.read.matrix
relative.abundance <-  final.read.matrix/rowSums(final.read.matrix)

final.otu.table <- otu_table(final.abundance, taxa_are_rows = FALSE)
relative.otu.table <- otu_table(relative.abundance, taxa_are_rows = FALSE)
tax.table <- tax_table(taxmat)

final.physeq <-  phyloseq(final.otu.table, tax.table)
relative.physeq <- phyloseq(relative.otu.table, tax.table)

save(final.physeq, file = snakemake@output[["physeq_abundance"]])
save(relative.physeq, file = snakemake@output[["physeq_rel_abundance"]])
write.table(sg_cluster, file = snakemake@output[["sg_cluster"]], row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
saveRDS(sg_reads, snakemake@output[["sg_reads"]])
