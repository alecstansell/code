
###################################################################################################################################################################################
#Plot sample tree to find obivous outliers
###################################################################################################################################################################################

remove_outliers_sample_tree <- function(clean_counts, sample_meta){
  
  #Transpose for function
  clean_counts <- t(clean_counts)
  #Run HAC clustering to plot tree
  sampleTree = hclust(dist(clean_counts), method = "average");
  
  #Select window size for plot
  sizeGrWindow(12,9)
  
  #Open pdf if you would like to write to pdf 
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  
  #Select appropriate params for axes and plot
  par(mar = c(0,4,2,0))
  par(cex = 0.6);
  
  #Plot clustered sample tree
  plot(sampleTree, main = "Sample Clustering ", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  
  #Decide whether to remove outliers
  print("Would you like to remove outliers? 
        1 for yes
        2 for no")
  outliers <- scan(file = "", nmax = 1)
  
  if(outliers == "1"){
    
    print("Select height to exclude sample outliers from")
    height <- scan(file = "", nmax = 1)
    
    # Plot a line to show the cut
    abline(h = height, col = "blue");
    
    print("View abline")
    invisible(readline(prompt="Press [enter] to continue"))
    
    # Determine cluster under the line
    clust = cutreeStatic(sampleTree, cutHeight = height, minSize = 10)
    
    #clust 1 contains the samples we want to keep.
    keepSamples = (clust==1)
    
    #Keep only samples that fall in the range we want
    datExpr = clean_counts[keepSamples, ]
    
    #Print number of samples
    print(paste("Samples are being excluded based on height of a height of", height, sep = " "))
    print(paste("You have", nrow(datExpr), "samples remaining", sep = " "))
    #atExpr = clean_counts
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    
    par(cex = 0.6);
    par(mar = c(0,4,2,0))
    
    sampleTree = hclust(dist(datExpr), method = "average");
    plot(sampleTree, main = "Sample Clustering ", sub="", xlab="", cex.lab = 1.5,
         cex.axis = 1.5, cex.main = 2)
    
    
    #Adjust meta data file to remove removed samples
    rownames(samples) <- samples$sample_id
    finalmeta <- samples[keepSamples,]
    
    input <- list(datExpr, finalmeta)
    
    return(input)
    
    
    
  } else{
    
    input <- ls(clean_counts, sample_meta)
    print("You chose to leave the counts unchanged")
    return(input)
    
    
  }
  
}


