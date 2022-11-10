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
colnames(taxonomy) <- c("Cluster","Taxonomy")


#setting important variables
gene_index <- seq(1,length(GeneLengths))
gene_names <- names(GeneLengths)
n.mapped.minimum <- 3 #The number of genes that needs reads that map to count the cluster as present
n.genes <- 100 # number of signature genes

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


# The read counts are normalized according to gene lengths
final.read.matrix <- matrix(NA, nrow=dim(Clusterlist[[1]])[2], ncol=length(Clusterlist))
sample.names <- colnames(Clusterlist[[1]]) 
rownames(final.read.matrix) <- sample.names
colnames(final.read.matrix) <- names(Clusterlist) 

final.Clusterlist <- Clusterlist

for (id in names(Clusterlist)){  
  # Repeat for the final SG
  final.gene.names <- screened_clusters[,3][screened_clusters[,'id']==id][1][[1]]$best
  final.colsum <- colSums(Clusterlist[[id]][final.gene.names, ])
  final.mapped <- names(final.colsum[final.colsum >= n.mapped.minimum])
  final.not_mapped <- sample.names[!sample.names %in% final.mapped]
  final.Clusterlist[[id]][,final.not_mapped] <- 0 # setting the counts to 0 if less than n.mapped.minimum genes have reads that map
  
  # The readcounts are divided by the gene length
  final.reads <- final.Clusterlist[[id]][final.gene.names, ] / GeneLengths[final.gene.names]
  
  # summing the read counts for the id/cluster/MGS
  final.read.matrix[, id] <- colSums(final.reads)
}

final.abundance <- final.read.matrix/rowSums(final.read.matrix)

final.otu.table <- otu_table(final.abundance, taxa_are_rows = FALSE)
tax.table <- tax_table(taxmat)

final.physeq <-  phyloseq(final.otu.table, tax.table)

save(final.physeq, file = snakemake@output[["physeq_abundance"]])
