set.seed(1337)

library(vcd) #goodfit
library(plyr)
library(VGAM) # for gpois
library(gridExtra) #for plotting ggplot next to each other
library(matrixStats) #rowMedians
library(parallel) #for parallizing; mclapply


# initializing

Clusterlist <- readRDS(snakemake@input[["clusters_sorted"]])
GeneLengths <- readRDS(snakemake@input[["gene_lengths"]])
Cluster <- readRDS(snakemake@input[["clusters_dir"]])

snakemake@source(snakemake@params[["functions"]])

#Number of signature genes
n.genes <- as.integer(snakemake@params[["n_genes"]])

#minimum number of mapped genes required
n.mapped.minimum <- as.integer(snakemake@params[["min_mapped_signature_genes"]])
n.minimum.samples <- as.integer(snakemake@params[["min_samples"]])

ids <- names(Clusterlist) #the ids of the MGS
id <- tail(strsplit(strsplit(snakemake@input[["clusters_dir"]], ".RDS")[[1]], "/")[[1]], n=1)

# initializing objects to keep track of stats throughout the process
mse.before <- rep(0, length(ids)) # MSE of the initial signature gene set according to the expected distribution
names(mse.before)  <- ids
mse.after <- (mse.before)  # MSE of the refined signature gene set
mse.final <- (mse.before)  # MSE of the refined signature gene set after removal of outliers

# identifying the names of genes in all MGS the dataset
#Cluster_gene_names <- c()
#for(cluster in names(Clusterlist)){
#  Cluster_gene_names <- c(Cluster_gene_names, rownames(Clusterlist[[Cluster]]))
#}

present_genes <- GeneLengths[names(GeneLengths) %in% rownames(Cluster)]


best_threshold <- c()
runtimes <- c()
mapped_samples <- c()
n_samples <- c()

run_one_id <- function(id){
  outputlist <- list()
  outputlist[['id']] <- id
  outputlist[['MSE']] <- list()
  outputlist[['genes']] <- list()

  #check that more than n genes pepresent in the id
  if (length(Clusterlist[[id]][,1])<=n.genes){
#    print(c("Not enough genes in the cluster. The number is", length(Clusterlist[[id]][,1])))
    return()
  }

  #colsum is the amount mapped genes in any given sample
  colsum <- colSums(Clusterlist[[id]][1:n.genes, ])
  # extracting only the samples which contains above n.mapped.minimum mapped reads -- if there are above n.mapped.minimum reads we believe this is a true detection
  mapped <- colsum[colsum >= n.mapped.minimum]
  
  n_samples <- setNames(c(n_samples, sum(colsum>=n.mapped.minimum)), c(names(n_samples), id))
  #if less than n.minimum.samples are mapped to the MGS, the MGS should be skipped
  if (length(mapped) < n.minimum.samples){
#    print("Not enough samples are mapped to the cluster")
    return()}
  genes <- names(Clusterlist[[id]][1:n.genes, 1])
  
  best.genes.step.zero <- genes
  ####################################################### The mean gene refinement #########################################################
  # loop over different threshold values for mean rank
  best.genes <- c()
  mse <- c()
  
  #fit without cutting anything out at all
  init.fit <- rank.mat(id = id,
                       gene.names = genes,
                       step = "mean",
                       phase = "initial",
                       threshold = 999,
                       t.start = 0,
                       t.end = n.genes,
                       mapped = mapped)
  
  #special case:
  #cant fix perfect. just stop here if you come to that 
  if(init.fit[['mse']]==0){
    outputlist[['MSE']][['initial']] <- 0
    outputlist[['MSE']][['mean']] <- 0
    outputlist[['MSE']][['best']] <- 0
    outputlist[['genes']][['initial']] <- init.fit[['good.genes']]
    outputlist[['genes']][['mean']] <- init.fit[['good.genes']]
    outputlist[['genes']][['best']] <- init.fit[['good.genes']]
    
    return(outputlist)
  }
  
  
  outputlist[['genes']][['init']] <- init.fit[["good.genes"]]
  outputlist[['mse']][['init']] <- init.fit[["mse"]]
  
  
  
  best.model <- c()
  best.model$mse <- init.fit[["mse"]]
  outputlist[["MSE"]][['initial']] <- init.fit[["mse"]]
  #a Finding the combination of genes which minimizes the MSE
  
  
  thresholds <- seq(35,60)
  t_best <- NA
  
  for (t in thresholds){
    i <- 0 #rounds of iterations
    #resets for new threshold
    genes <- names(Clusterlist[[id]][1:n.genes, 1])
    MSE.old <- best.model[['mse']]
    MSE <- ""
    new_genes_assigned <- F
    
    #carries out the initial fit
    init.genes <- rank.mat(id = id,
                           gene.names = genes,
                           step = "mean",
                           phase = "initial",
                           threshold = t,
                           t.start = 0,
                           t.end = n.genes,
                           mapped = mapped)
    
    # iteratively replace genes until the MSE does not decrease (with an interval for MSE allowed)
    # if the mse is seen before, break
    
    while ((MSE < (MSE.old * 1.01)) & (! MSE %in% mse[1:length(mse) - 1])){ # if the mse is seen before, break
      
      #If there has been new genes assigned, then you can re-perform an initial fit
      if(new_genes_assigned==T){
        init.genes <- rank.mat(id,
                               gene.names = genes,
                               "mean",
                               "initial",
                               t,
                               0,
                               n.genes,mapped = mapped)
      }
      
      if (length(init.genes$good.genes) < 10){ # if there is below 10 good genes it is impossible to fit model and there is no reason to keep going, hence the break
        MSE <- init.genes$mse
        break
      }

      if ((min(500, length(Clusterlist[[id]][, 1]))-n.genes)<(n.genes-length(init.genes$good.genes))){ ##check if there is enough potential replacement genes in the geneset
        MSE <- init.genes$mse
        best.model <- init.genes
        best.genes <- genes
        break
      }
      
      #if all the genes meet the threshold set up before then we are done
      #if we have 100 good genes then there is no reason to keep going
      if (length(init.genes$good.genes) == n.genes){
        MSE <- init.genes$mse
        best.model <- init.genes
        best.genes <- init.genes$good.genes
        
        break
      }
      
      "if this is the first iteration, then set the MSE.old as the current MSE"
      if (i == 0) {
        MSE.old <- init.genes$mse
        mse.before[id] <- init.genes$mse
        
        #if not, then update the old MSE as the MSE of the from the previous iteration
      } else {
        MSE.old <- rotation$mse
      }
    if ((n.genes+1)==min(500, length(Clusterlist[[id]][, 1]))){ # if there is only one gene to rotate it is not enough
        MSE <- init.genes$mse
        best.model <- init.genes
        best.genes <- init.genes$good.genes
        break
      } 
    
      rotation <- rank.mat(id = id,
                           gene.names = init.genes$good.genes,
                           step = "mean",
                           phase = "rotation",
                           threshold = t,
                           t.start = n.genes,
                           t.end = min(500, length(Clusterlist[[id]][, 1])),
                           n.replace = (n.genes - length(init.genes$good.genes)),
                           mapped = mapped)
      
      
      # if the initial gene set is better than the rotation, then why continue? <You started out with the best thing
      if (init.genes$mse < rotation$mse & i == 0){
        MSE <- init.genes$mse
        best.model <- init.genes
        best.genes <- genes
        break
      }
      MSE <- rotation$mse
      rotation[['threshold']] <- t
      
      # best.genes <- c(init.genes$good.genes, names(rotation$good.genes))
      i <- i + 1
      mse <- c(mse, MSE.old)
      
      # Keeping the best model across iterations and thresholds if the new model is better than the best model across thresholds
      if ((MSE <= MSE.old) & (MSE < best.model$mse)){
        tmp.best.genes <- c(init.genes$good.genes, names(rotation$good.genes))
        if (length(tmp.best.genes)<100){
          print("There is not enough genes to rotate")
          break
        }
        best.genes <- tmp.best.genes  
        best.model <- rotation
        genes <- best.genes
        t_best <- t

        #if you assign new genes, you can do another initial fit
        new_genes_assigned <- T
      }else{
        new_genes_assigned <- F
      }
    }#while over
  }
  
  outputlist[["genes"]][['mean']] <- best.genes
  outputlist[['MSE']][['mean']] <- best.model[['mse']]
  
  best.genes.step.one <- best.genes
  
  
  ############################################################## The 95-percentile gene refinement ##########################################################3
  #This step does: checks the best fits are still okay
  best.old <- best.genes
  best.old.mse <- best.model$mse
  thres <- c(90,91,92,93,94,95,96,97,98) # The thresholds for 95-percentile ranks
  
  #if you never found a new better set, then you can just start over
  if(is.null(best.old)){
    best.old <- rank.mat(id,
                         gene.names = genes,
                         "mean",
                         "initial",
                         999,
                         0,
                         n.genes,mapped = mapped)[['good.genes']]
  }
  
  t_best <- NA
  for (t in thres){
    MSE <- 0
    MSE.old <- best.old.mse
    k <- 0 #k is just the new i (iteration counter)
    new.genes <- best.old
    
    #initially, evaulated the init.quant
    genes_revised <- T 
    
    while((MSE < (MSE.old * 1.05)) & (! MSE %in% mse[1:length(mse) - 1])){
      
      if(genes_revised==T){
        init.quant <- rank.mat(id,new.genes, "quantile", "initial", t, 0, n.genes,mapped=mapped)
      }
      
      if (k != 0){
        MSE.old <- rot.quant$mse
      }else {
        MSE.old <- init.quant$mse
      }
      #redoes the fit with new genes found during rotation, if there had been found new during rotation
      if(genes_revised==T){
        init.quant <- rank.mat(id,new.genes, "quantile", "initial", t, 0, n.genes,mapped=mapped)
      }
      
      # If the initial genes are all included by the threshold, then we are good and break
      if (length(init.quant$good.genes) == n.genes){
        MSE <- init.quant$mse
        if (MSE < best.model$mse){
          best.model <- init.quant
          best.genes <- init.quant$good.genes}
        break
      }
      
      if ((min(500, length(Clusterlist[[id]][, 1]))-n.genes)<(n.genes-length(init.quant$good.genes))){ # if there is not more new genes than genes required for rotation
        best.model <- init.quant
        best.genes <- new.genes
        break
      }

    if ((n.genes+1)==min(500, length(Clusterlist[[id]][, 1]))){ # if there is only one gene to rotate it is not enough
        MSE <- init.genes$mse
        best.model <- init.genes
        best.genes <- init.genes$good.genes
        break
      }

      # if there is only 5 good genes it is impossible to fit model
      if (length(init.quant$good.genes) < 5){ # if there is only 5 good gene it is impossible to fit model
        if (init.quant$mse < best.model$mse){#ØP: and if that model, using 5 or fewer, is still better, then use that?
          best.model <- init.quant
          best.genes <- new.genes
        }
        #print(c("There are not enough good genes to fit the model, the numer of good genes is:", length(init.quant$good.genes)))
        break
      }
      
      #rotates
      rot.quant <- rank.mat(id,
                            init.quant$good.genes,
                            "quantile",
                            "rotation",
                            t,
                            n.genes,
                            min(500, length(Clusterlist[[id]][, 1])),
                            (n.genes - length(init.quant$good.genes)),
                            mapped = mapped)
      new.genes <- c(init.quant$good.genes, names(rot.quant$good.genes))
      MSE <- rot.quant$mse
      k <- k + 1
      mse <- c(mse, MSE)
      
      if (MSE < best.model$mse){
        tmp.best.genes <- c(init.quant$good.genes, names(rot.quant$good.genes))
        if (length(tmp.best.genes)<100){
          print("There is not enough genes to rotate")
          break
        }
        best.genes <- c(init.quant$good.genes, names(rot.quant$good.genes))
        best.model <- rot.quant

        genes_revised <- T
        t_best <- t
      }else{
        genes_revised <- F
      }
    }
  }
  
  outputlist[["genes"]][['best']] <- best.genes
  outputlist[["MSE"]][['best']] <- best.model$mse
  
  return(outputlist)
}


output <- run_one_id(id)

saveRDS(output, file=snakemake@output[["cluster_screened"]])
