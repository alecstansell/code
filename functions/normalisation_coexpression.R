###########################################################################################################################################################################################
# Normalissation Code based off normalisation_counts.R (by Alec)
###########################################################################################################################################################################################

#Load the DESeq2 package
suppressPackageStartupMessages(library("DESeq2"))
# input_counts <- raw_counts
# input_meta <- samples
# 


normalise_it <- function(input_counts, input_meta){
  
  #Meta data cannot have rownames for DESeq
  rownames(input_meta) <- NULL
  
  #Create object with all gene names
  genes <- rownames(input_counts)
  
  #Convert geneCounts data frame to matrix
  geneCounts <- as.matrix(input_counts)
  
  #Convert counts in geneCounts matrix to integers by column
  countMatrix <- (apply(geneCounts, MARGIN = 2, FUN = as.integer))
  
  rownames(countMatrix) <- genes
  
  
  #Genes have to be unique for use in Deseq2 check if they are
  if(all(isUnique(genes)) == FALSE){
    stop("Gene names are not unique. You can add unique gene identifies using the make.unique(genes) command")
  }
  
  #Create dds in order to run normalisation
  dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = input_meta, design = ~ 1)
  
  #Run deseq2 note that there is no design condition as we are only looking to normalise the data
  dds <- DESeq(dds)
  
  #Get results
  res <- results(dds)
  
  #Viewsize factors - ideally they should be close to 1
  sizeFactors(dds)
  
  #Create object containing the size factors
  sizefactor_counts <- sizeFactors(dds)
  
  #Calculate colsums to compare between samples - large differences indicate differences in depth
  colSums(counts(dds))
  
  #Plot column sums according to size factor
  plot(sizeFactors(dds), colSums(counts(dds)))
  abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
  
  
  #The argument normalized equals true, divides each column by its size factor.
  size_factor_counts <- log2( counts(dds, normalized=TRUE) + 1 )
  
    #Variance stabilised counts
  variance_stabilised_counts <- getVarianceStabilizedData(dds)
  
  ## R log transformation
  rlog_counts <- rlogTransformation(dds)
  rlog_counts <- assay(rlog_counts)
  
  #Add names a back to rlog counts
  rownames(rlog_counts) <- genes
  
  #Log transformationt to make data more "microarray like"
  log2_input_counts <- log2(input_counts + 1)
  
  #Create large normalised list with the 4 different options
  all_counts <- list(size_factor_counts, variance_stabilised_counts, rlog_counts, input_counts, input_meta)
  names(all_counts) <- c("size_factor_counts", "variance_stabilised_counts", "rlog_counts", "input_counts", "Input Metadata")
  
  return(all_counts)
  
  
  
}
  



