library(ggplot2)
library(vcd) #goodfit
library(plyr)
library(VGAM) # for gpois
library(matrixStats) #rowMedians



#dummy variables for testing purposes

if(1==0){
  Genes <- abs(rnorm(100,1,1))
  names(Genes) <- as.character(seq(1,100))
  Reads <- Genes+1
  mse=10
  
  
  gene.prob <- runif(n = 100,min = 0,max = 1)
  names(gene.prob) <- seq(1,100)
  gene.pos <- seq(1,100)
  names(gene.pos) <- seq(1,100)
  
  
  
  gene.names <- init.genes$good.genes
  step='mean'
  phase='rotation'
  threshold=t
  n.genes=n.genes
  t.end=min(500, length(utuseq_MGSgeneCountsML[[id]][, 1]))
  t.start=(n.genes - length(init.genes$good.genes))
  mapped <- mapped
  
}

# loading functions
pois.plot <- function(Genes, Reads, mse, id){
  # Plotting the reads and mapped genes for each sample
  #
  # Args:
  #   Genes: A named vector containing the number of genes mapped in the samples
  #   Reads: Vector containing the number of reads mapping to the Genes. Genes and Reads must 
  #           have the same length, greater than one, with no missing values
  #   mse:  The MSE of the MGS when compared to the expected distribution
  #   id is the ID of the species in the HGMGS object
  #
  # Returns:
  #   A plot of all samples with their corresponding genes and mapped reads (log)
  
  
  df_plot <- do.call(rbind, Map(data.frame, Reads = Reads, Genes = Genes, name = names(Genes)))
  p = ggplot(log = "x", xlim = c(0.1, 1e5), cex = 0.8, pch = 25, ylim = c(0, n.genes), 
             data = df_plot, aes(text = name, y = Genes, x = Reads)) + 
    geom_point(size = 0.7, aes(color = "Samples")) +
    stat_function(fun =  function(Genes) (1 - ((n.genes - 1) / n.genes)^Genes) * n.genes, aes(Genes, color="Predicted"), size=0.2, alpha=0.5) + #t he expected distribution
    xlab('Reads mapped (log)') +
    ylab('Signature genes detected') +
    ggtitle(id)
  
  plot.data <- (p + scale_x_log10() + scale_colour_manual(values = c("Samples" = cmcol1[2], "Outlier" = cmcol1[5], "Predicted"="black")))
  plot.data <- plot.data + annotate("text", x=20, y=80, label = paste("MSE:", round(mse, 1)), parse=T, size=3)
  plot.data <- plot.data  + theme_minimal() + theme(legend.position="bottom")
  plot.data  # for creating PDF
} 


gene.selection <-  function(n.replace, gene.prob, gene.pos) {
  #  Identifying the best suited replacement signature genes
  # Args:
  #   n.replace: number of genes to replace
  #   gene.prob: Vector with the probability and name of the genes
  #   gene.pos: Vector with the position of the gene
  # Returns:
  #   replace: Vector with genes with a high Pearson residual in one or more samples. 
  
  if (n.replace > length(gene.prob)){
    n.replace <- length(gene.prob)
  }
  
  replace.genes <- sort(abs(gene.prob) * 0.975 + gene.pos * 0.025, decreasing = FALSE)[0:n.replace]
  
  return("replace" = replace.genes)
}

rank.mat <- function(id, gene.names, step, phase, threshold, t.start, t.end, n.replace = FALSE, mapped) {
  #  Identify suitable signature genes by fitting NB model and ranking the genes
  #
  # Args:
  #   gene.names: A list of genes used to train the model
  #   step: Indicating whether the suitable genes are found based on "median" or "quantile"
  #   phase: The phase of the modelling. Currently there are two phases: 
  #     "initial": fit NB model to the original signature genes
  #     "rotation": fit NB model to refined signature genes 
  #   threshold: Mean value of rank for keeping genes
  #   t.start: The start of the genes in the read count matrix, default: 0 for initial, 100 for rotation
  #   t.end: The end of the genes in the read count matrix, default: 100 for initial, 500 for rotation
  #   n.replacee: The number of genes to replace, int
  #   mapped: samples we are confident the MGS is present in (e.g. have more than n.mapped.minimum reads mapping)
  #
  # Returns:
  #   rank.mat: Matrix containing the ranks of the gene.names in the samples
  #   good.genes: List of the gene names of the genes that we keep in the signature gene set, for rotation it is the list of genes of length n.replace
  #   countReads: double with sample name and corresponding total read count
  #   countGenes: double with sample name and corresponding detected number of signature genes
  #   mse: The MSE of the fit to the expected distribution (s=(1-((G-1)/G)^N)*G)
  
  rank.matrix <- matrix(NA, nrow=(t.end-t.start), ncol=length(mapped))
  
  reads <- round(Clusterlist[[id]][gene.names, names(mapped)] / 
                   (present_genes[rownames(Clusterlist[[id]][gene.names, ])] * 10^-3))
  
  # running goodfit on all samples, saving the model in the gf matrix
  gf.nbin <- apply(matrix(1:length(mapped), nrow = length(mapped), ncol = 1), 1, function(x) (goodfit(reads[, x], type = "nbinom", method = "ML", par = list(size = mean(reads[,x])))))
  
  
  # filling in the rank.matrix with the rank of all genes within each sample
  if (phase == "initial") {
    rank.matrix <- mapply(function(c) rank.matrix[, c] = rank((resid(gf.nbin[[c]])[gf.nbin[[c]]$count[reads[,c]+1]+1]), ties.method = "average"), 1:length(mapped)) 
    
  } else if (phase == "rotation") {
    new.reads <- matrix(NA, nrow = (t.end-t.start), ncol = length(mapped))
    new.reads <- round(Clusterlist[[id]][(t.start + 1):t.end, 1:length(names(mapped))] / 
                         (present_genes[rownames(Clusterlist[[id]][(t.start + 1):t.end, ])] * 10^-3))
    
    # if a new gene have more reads than found in the good genes, then set the readcount to the max observed in the good genes
    for (i in 1:length(mapped)) new.reads[, i][new.reads[, i] > max(gf.nbin[[i]]$count)] = max(gf.nbin[[i]]$count)
    rank.matrix <- mapply(function(c) rank.matrix[, c] = rank((resid(gf.nbin[[c]])[gf.nbin[[c]]$count[new.reads[, c]+1]+1]), ties.method = "average"), 1:length(mapped)) 
  }
  
  rownames(rank.matrix) <- names(Clusterlist[[id]][(t.start + 1):t.end, 1]) # genes as rownames
  colnames(rank.matrix) <- names(mapped) # samples as colnames
  
  #######################################################################
  # Adding step to weigh samples with more reads higher
  ######################################################################
  ori_mean = mean(rowMeans(rank.matrix)) #50.5 
  min = min(rowMeans(rank.matrix))- 5
  max = max(rowMeans(rank.matrix)) + 5
 
  rank.matrix <- t(t(rank.matrix) * log10(colSums(reads))) 
  rank.matrix <- replace(rank.matrix, rank.matrix == -Inf, 0)
  
  scale = ori_mean/mean(rowMeans(rank.matrix))
  rank.matrix <- replace(rank.matrix*scale, rank.matrix*scale == 'NaN', 0)

  
  if (phase == "initial") {
    
    # identifying the genes to replace
    if (step == "mean"){
      gene.performance <- rowMeans(rank.matrix) # identifying the means of the signature genes
      
    } else if(step=="quantile"){
      gene.performance <- rowQuantiles(rank.matrix, probs = seq(from = 0, to = 1, by = 0.05))[,20] #identifying the 95-percentile of the signature genes
    }
    
    low.prob.genes <- gene.names[gene.performance > threshold] # if the rank of the genes are above threshold
    n.replace <- length(low.prob.genes) # number of genes to replace
    good.genes <- gene.names[!(gene.names %in% low.prob.genes)] # the genes with rank below the threshold
    
    genes_r <- Clusterlist[[id]][c(names(good.genes),gene.names), ]
    
    
    final.reads <- round(genes_r / (present_genes[rownames(genes_r)] * 10^-3))

    countGenes <- colSums(final.reads > 0) 
    countReads <- round(colSums(final.reads))
    
  } else if (phase == "rotation"){
    if (step == "mean"){
      gene.performance <- rowMeans(rank.matrix)
    } else if(step=="quantile"){
      gene.performance <- rowQuantiles(rank.matrix, probs = seq(from = 0, to = 1, by = 0.05))[,20] #95percentile
    }
    
    names(gene.performance) = names(Clusterlist[[id]][(t.start + 1):t.end, 1]) 
    high.rank <- gene.performance 
    high.pos <- seq((t.start + 1), t.end)
    
    # ensuring that the new genes are not already part of the signature gene set
    high.rank <- high.rank[! names(gene.performance) %in% gene.names]
    high.pos <- high.pos[! names(gene.performance) %in% gene.names]
    
    # identifying the best replacement genes
    good.genes <- gene.selection(n.replace, high.rank, high.pos)
    
    genes_r <- Clusterlist[[id]][c(names(good.genes),gene.names), ]
    
    final.reads <- round(genes_r / (present_genes[rownames(genes_r)] * 10^-3))

    countGenes <- colSums(final.reads > 0) 
    countReads <- round(colSums(final.reads))
    
  }
  
  # calculating the predicted gene counts
  pred <- (1 - ((n.genes - 1) / n.genes)^countReads) * n.genes
  
  # calculating the MSE
  mse <- mean((countGenes - pred)^2)
  
  return(list("good.genes" = good.genes, "countReads" = countReads, "countGenes" = countGenes, "mse" = mse,'countReads'=countReads,'countGenes'=countGenes))
}

