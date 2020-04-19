#Load limma
library(limma)
library(edgeR)
source('~/doorknobdave/alecS/alec_code/functions/remove_duplicates.R')
source('~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/distribution_plot.R')

unlock_the_limma <- function(raw_counts, meta_data){
  
 #Prepare counts
  d0 <- DGEList(raw_counts)
  d0 <- calcNormFactors(d0)
  cutoff <- 6
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 
  dim(d)
  
  # Remove genes with zero variance across samples - this helps the model
  #limma_input <- raw_counts[apply(raw_counts, 1, var) > 0,]
  
  #Plot sample distribution
  distributionplot(d, "removed") + xlim(0, 10)
  
  invisible(readline(prompt="Press [enter] to continue"))
  
  #### Use only differentially expressed genes for coexpression analysis to yield the most informative results
  #Batch effects can be accounted for here - my samples dont seem to contain any batch clustering upon analysis with sleuth PCAs
  
  #Example of this in the temr of linear model of limma
  # mod <- model.matrix(~0+samples$cell+samples$batch)
  
  #Design matrix for Limma differential expression analysis 
  mod <- model.matrix(~0+meta_data$condition)
  
  file_out_limma <- "~/doorknobdave/alecS/alec_code/degs_output/healthy_cells/limma"
  
  invisible(readline(prompt="Press [enter] to continue"))
  
  
  #Voom transformation and make mean variance plot
  png(filename=file.path(file_out_limma, "mean_variance_plot.png"), width = 600, height = 400, )
  voomtrans <- voom(d, mod, plot = T)
  dev.off()
  
  invisible(readline(prompt="Press [enter] to continue"))

  #Change column names of mod to make it easier to understand
  colnames(mod) <- levels(meta_data$condition)
  
  #Fit modelin limma
  fit <- lmFit(voomtrans, design=mod)
  
  #Generate a list of all possible pairwise contrasts in my case this is only 3
  cell_pairs <- t(combn(levels(meta_data$condition), 2))          
  
  #Create a list of all the comparisons  
  comparisons <- list()               
  for (i in 1:nrow(cell_pairs)) {                                                                                                                                     
    comparisons[[i]] <- as.character(cell_pairs[i,])                                                                                                      
  }    
  
  #Create vector to store sig genes
  sig_genes <- c()
  
  #Print possible compariosns 
  print(comparisons)
  
  #Which comparison would you like
  print('Pick Comparison 1, 2 or 3')
  comp <- scan(file = "", nmax = 1)
  
  
  #Set contrast formula
  contrast_formula <- paste(comparisons[[comp]], collapse=' - ')
  
  #Generate contrast Matrix
  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=mod)
  
  #Fit contrasts
  contrast_fit <- contrasts.fit(fit, contrast_mat)
  
  #Run
  eb <- eBayes(contrast_fit)

  #Extract Sig genes for comparison
  #topTable(eb, number=Inf, p.value=0.05)
  
  #Retrieve sig genes
  genes <- topTable(eb, number=Inf)
  genes <- as.data.frame(genes)
  sig_genes <- genes[which(genes$P.Value < 0.05),]
  
  #Get names
  sig_genes$Gene_Name <- rownames(sig_genes) 
  
  #View length of the genes
  print(paste("The number of genes is", length(rownames(sig_genes)), sep = " "))
  
    #Remove duplicates with remove dup function
  sig_genes <- removedup(sig_genes)
  rownames(sig_genes) <- sig_genes$Gene_Name
  
  print(paste("The number of genes after duplicates are removed is", length(rownames(sig_genes)), sep = " "))
  
  #return sig genes for comparison
  return(sig_genes)
  
}
