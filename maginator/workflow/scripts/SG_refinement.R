set.seed(1337)
window_size <- 699
turkey_factor <- 1.2

#loading CM colors
cmcol1 = c('#0571B0', '#92C5DE', '#999999', '#E0E0E0', '#CA0020',
           '#FF8784', '#FFD900', '#FFE967', '#349F2C', '#BAE4B3',
           '#6B3D9A', '#CAB2D6')

###############################################################################
##################################### MAIN ####################################
###############################################################################
library(ggplot2)
library(vcd) #goodfit
library(plyr)
library(VGAM) # for gpois
library(gridExtra) #for plotting ggplot next to each other
library(matrixStats) #rowMedians
library(parallel) #for parallizing; mclapply


# initializing

Clusterlist <- readRDS(snakemake@input[["R_clusters"]])
GeneLengths <- readRDS(snakemake@input[["R_gene_lengths"]])

snakemake@source(snakemake@params[["functions"]])

#Number of signature genes
n.genes <- 100

#minimum number of mapped genes required
n.mapped.minimum <- 3

ids <- names(Clusterlist) #the ids of the MGS

# initializing objects to keep track of stats throughout the process
mse.before <- rep(0, length(ids)) # MSE of the initial signature gene set according to the expected distribution
names(mse.before)  <- ids
mse.after <- (mse.before)  # MSE of the refined signature gene set
mse.final <- (mse.before)  # MSE of the refined signature gene set after removal of outliers

# identifying the names of genes in all MGS the dataset
Cluster_gene_names <- c()
for(Cluster in names(Clusterlist)){
  Cluster_gene_names <- c(Cluster_gene_names, rownames(Clusterlist[[Cluster]]))
}

present_genes <- GeneLengths[names(GeneLengths) %in% Cluster_gene_names]

# looping over the ids of the MGS
#registers time to run each
best_threshold <- c()
runtimes <- c()
mapped_samples <- c()
n_samples <- c()


if(file.exists(snakemake@output[["screened_clusters"]])){
  load(snakemake@output[["screened_clusters"]])
   Clusterlist <- clusterlist_screened
   rm(clusterlist_screened)
}else{
  #preselecting genes
  t0 <- Sys.time()
  print('pre-screeing genes')
  for(id in ids){

    genes_r <- Clusterlist[[id]]
    final.reads <- round(genes_r / (present_genes[rownames(genes_r)] * 10^-3))

    #finds the frequency of each gene. We want a ´p_sig_gene~1/n_signature_gene
    #Because if we just have one that's all messy and much more likely than the others, then you effed. having a high one is actually a lot worse, since 95 or more is cool
    frequency_of_genes <- rowSums(final.reads>0, na.rm = T)

    #So! Basically...
    #if I doon't exclude certain genes, I end up with scenarios in which U get a 3-4 genes that are 5x as likely to be identified compared to the ones that are hard to find
    #That means, that a LOT of instances where if you sequence, say, 10 times, anything below nine is significant given a totally random draw. But!
    #If you have 3-4 genes that are very likely to pop up, you end up getting a ton of times where the same high-frequency gene gets found TWICE!
    #This is exemplified in MGS.hg0021, where, almost all of the high-read items gets has 98/100, indicating that there are probably no instances of
    #like, species that lack a substantial part of their signature genes, so nothing thrown away on the far right end of the plot
    #however, in the "climbing" part of the plot, a ton of samples are labeled as bad. If you go in and look at which genes are found how many times here, you often find that
    #a few genes soak up a lot of the reads, leaving the others lacking :(
    #So, we want to choose a set of genes we can purify on that 1) are usually found in high frequency and 2) are relatively consisent with each other.

    #in a test, I tried to select all genes between median*1.2 and median/1.2 and only run on those. Bad sutff. a lot of high-read samples never got all genes, good stuff: they followed a n_expected_signature_genes =80 extremely well
    #But, of course, I did choose a lot of genes that are sort of middling in terms of frequency, i.e. relatively hard to detect, indicating that they probably aren't found in all members of that MGS...
    #So, let's try to only select high-frequency genes, none of which are outliers by turkey's method, but a bit more sensitive, because I think I get get away with it basically...
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
    i <- match(id,ids)
    print(paste0(round(i/length(ids)*100),'% screened. Approximate time left:'))

    time_taken_per_sample <- (Sys.time()-t0)/i
    time_left <- length(ids)-i*time_taken_per_sample
    print(time_left)
  }
  clusterlist_screened <- Clusterlist
  save(clusterlist_screened, file = snakemake@output[["screened_clusters"]])
}


utuseq_MGSgeneCountsML <- Clusterlist





n_id=0 
t0=Sys.time()


run_one_id <- function(id){
  outputlist <- list()
  outputlist[['id']] <- id
  outputlist[['MSE']] <- list()
  outputlist[['genes']] <- list()
  message(id)

  #check that more than n genes pepresent in the id
  if (length(Clusterlist[[id]][,1])<=n.genes){
    print(c("Not enough genes in the cluster. The number is", length(Clusterlist[[id]][,1])))
    return()
  }

  #colsum is the amount mapped genes in any given sample
  colsum <- colSums(Clusterlist[[id]][1:n.genes, ])
  # extracting only the samples which contains above 3 mapped reads -- if there are 3 reads we believe this is a true detection
  # ØP: shouldn't it be something more like "how many genes we think we can detect reliably?" Surely there is a difference between detecting 4 genes once and detecting one gene four times?
  mapped <- colsum[colsum >= n.mapped.minimum]
  
  n_samples <- setNames(c(n_samples, sum(colsum>=n.mapped.minimum)), c(names(n_samples), id))
  #if less than 3 samples are mapped to the MGS, the MGS should be skipped
  if (length(mapped) < 3){
    print("Not enough samples are mapped to the cluster")
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
    # any reason for the 1% decrease requirement?
    
    #TRINE HAS CHANGED: while ((MSE < (MSE.old * 1.01)) & (! MSE %in% mse[1:length(mse) - 1])){ # if the mse is seen before, break
    
    while ((MSE < (MSE.old * 1.01)) & (! MSE %in% mse[1:length(mse) - 1])){ # if the mse is seen before, break
      #actually does the fit, test and stuff
      #ØP: Wait is there a reason for redoing this for every iteration? Can't we just recycle the old iteration if we don't update the geneset through an improved MSE? Seems like a lot of time is wasted here?
      
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
      
      #ØP: Is there a particular reason for you including countReads and countGenes? Maybe some time could be shaved off by not having to export all those vectors of length nsample...
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
  
  # print(c("This is the best mse after mean", best.model$mse))
  ############################################################## The 95-percentile gene refinement ##########################################################3
  #You really should look at how well stuff fits at its best! This is what this step does: checks the best fits are still okay
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
 #   print(t)
    MSE <- 0
    MSE.old <- best.old.mse
    #ØP: k is just the new i (iteration counter)
    k <- 0
    new.genes <- best.old
    
    #initially, evaulated the init.quant
    genes_revised <- T #F
    
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
      
      #ØP: if the initial genes are all included by the threshold, then we are good and break
      if (length(init.quant$good.genes) == n.genes){
        #print("The initial genes are all good")
        MSE <- init.quant$mse
        if (MSE < best.model$mse){
          best.model <- init.quant
          best.genes <- init.quant$good.genes}
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
        best.genes <- c(init.quant$good.genes, names(rot.quant$good.genes))
        if (length(best.genes)<100){
          print("There is not enough genes to rotate")
          break
        }
        
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
  
  
  i <- match(id,ids)
  time_per_sample <- (Sys.time()-t0)/i
  samples_remaining=length(ids)-i 
  
  return(outputlist)
}


# The signature gene optimization is run
numCores <- snakemake@threads #the number of threads used are extracted from the Snakefile

outputs_parallel <- mclapply(ids, run_one_id, mc.cores = numCores)

names(outputs_parallel) <- ids

saveRDS(outputs_parallel, file=snakemake@output[["screened_clusters"]])

gene_index <- seq(1,length(GeneLengths))
names(gene_index) <- names(GeneLengths)

gene_names <- names(gene_index)
names(gene_names) <- unname(gene_index)


trine_MGS_object <- list()
trine_MGS_object[['i']] <- list()

allMGS <- names(outputs_parallel)
nMGS <- length(allMGS)

i <- 0

for(MGS in allMGS){
  i=i+1

  all_genes <- which(names(GeneLengths) %in% rownames(Clusterlist[[MGS]])) # getting all gene index numbers of the Cluster

  
  for(entry in names(outputs_parallel[[MGS]][["genes"]])){
    sig_genes <- which(names(GeneLengths) %in% outputs_parallel[[MGS]][["genes"]][[entry]]) 
    trine_MGS_object[['i']][[paste(MGS,entry,sep='_')]] <- c(sig_genes, all_genes[!all_genes%in%sig_genes]) 

  }
  if(length(all_genes)>100){
  random_sig_genes <- sample(all_genes, 100)}
  else{random_sig_genes <- sample(all_genes, length(all_genes))}
  trine_MGS_object[['i']][[paste(MGS,'random',sep='_')]] <-c(random_sig_genes, all_genes[!all_genes%in%random_sig_genes])
}



save(trine_MGS_object, file=snakemake@output[["MGS_object"]])
