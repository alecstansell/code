##############################################################################################################################################################################################
#Master Self Organising Map Script
##############################################################################################################################################################################################

##############################################################################################################################################################################################
#Adapted and written by Alec Stansell
#Shane Lynn 14-01-2014 
##############################################################################################################################################################################################

#Make sure results are reproducible
set.seed(1)

##############################################################################################################################################################################################
#Load Packages
##############################################################################################################################################################################################

#List packages
packages <- c("kohonen", "dummies", "ggplot2", "limma", "sp", "reshape2", "maptools")

#install if necessary
#install.packages(packages)

#Load packages
lapply(packages, require, character.only = TRUE)

##############################################################################################################################################################################################
#Input datasets
##############################################################################################################################################################################################

# Colour palette definition
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')
mypalette <-c("#345eeb", '#7134eb', '#3bcc40', '#3bb1cc', '#cc423b','#cc423b')

#Set working directory
setwd('~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/data')

#Load normalised count data for use in SOM
load("final_heat_normalised.RData")

#Set specific normalisation version
counts_input <- normalised_clean$variance_stabilised_counts


##############################################################################################################################################################################################
#SOM return cluster function
##############################################################################################################################################################################################

som_return_cluster <- function(user_interface = FALSE, cluster, rlen){
  
  #Get file names for DESeq2 differentially expressed genes
  file_name <- list.files("~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/pairwise_no_dups")
  
  #Choose dataset to use
  
  print("Choose file 1, 2, or 3")
  print(file_name)
  
  choice <- scan(file = "", n = 1)
  
  file_choice <- file_name[choice]
  
  # Import DE genes from file if you are not going to use the above
  sig <- read.csv(file.path("~/doorknobdave/alecS/alec_code/coexpression/final_coexpression/pairwise_no_dups", file_choice))
  sig <- sig[,-1]
  colnames(sig) <- c("gene_names", "log2fold")
  
  #Set which cells to be used (so as to subset large data set for comparison)
  
  if(file_choice == "EMvsCM_dup_rem_sig.csv"){
    
    cell1 <- "EM"
    cell2 <- "CM"
    
  } else if(file_choice == "NvsCM_dup_rem_sig.csv"){
    
    cell1 <- "N"
    cell2 <- "CM"
    
    
  } else if(file_choice == "NvsEM_dup_rem_sig.csv"){
    
    cell1 <- "N"
    cell2 <- "EM"
    
  }
  
  #Pick cells for comparison
  som_input <- counts_input[, which(grepl(cell1, colnames(counts_input)) | grepl(cell2, colnames(counts_input)) == TRUE)]
  
  #Choose whether you would like to select specific upregulated or downregulated genes for coexpression analysis
  if(user_interface == TRUE){
    
    #Print choices
    print("Choose whether you would like to use downregulated, upregulated or all DE genes")
    
    #Press enter
    invisible(readline(prompt="Press [enter] to continue"))
    
    #Print selection
    print("Type up for upregulated, down for downregulated and all for all Differentially expressed genes")
    
    #Read selection
    reg <- scan(file = "", what = "character", nmax = 1)
    
  }else{
    
    reg <- 'all'
    
  }
  
  #If statement to select specific genes needed 
  if(reg == "up") {
    df <- subset(sig, log2fold > 0)
    df <- df[order(df$log2fold, decreasing = FALSE),]
    print("Upregulated genes selected")
  } else if(reg == "down"){
    df <- subset(sig, log2fold < 0)
    df <- df[order(df$log2fold, decreasing = TRUE),]
    print("Downregulated genes selected")
  } else if(reg == "all"){
    df <- sig
    df <- df[order(df$log2fold, decreasing = TRUE),]
    print("All genes selected")
    
  }else{
    
    print("Are you fucking retarded")
  }
  
  #Subset counts for your specific comparison
  som_input <- som_input[rownames(counts_input) %in% df$gene_names,]
  
  ##############################################################################################################################################################################################
  #Clean data
  ##############################################################################################################################################################################################
  
  #remove incomplete samples:
  incompletes <- which(!complete.cases(som_input))
  length(incompletes)
  #where the avr_education_level is NaN - replace with mean
  # data$avr_education_level[incompletes] <- mean(data$avr_education_level, na.rm=TRUE)
  #recalculate after adjustment
  incompletes <- which(!complete.cases(som_input))
  if (length(incompletes) > 0){
    print(paste0("Removing ", length(incompletes), " data points that have missing values."))
    som_input <- som_input[-incompletes, ]
  }
  rm(incompletes)
  
  ##############################################################################################################################################################################################
  #Train Self Organising Map
  ##############################################################################################################################################################################################
  
  #choose the variables with which to train the SOM
  #Naive data train
  data_train <- som_input
  
  #Train SOM using Kohonen method
  data_train_matrix <- as.matrix(scale(data_train))
  
  #Set names
  names(data_train_matrix) <- names(data_train)
  
   print("Training SOM model")
  
  # SOM Model
  system.time(som_model <- som(data_train_matrix, 
                               grid=som_grid, 
                               rlen=rlen, 
                               alpha=c(0.1,0.01), 
                               keep.data = TRUE ))
  #Plot training progress
  print("Plotting training progress for quality evaluation")
  plot(som_model, type = "changes")
  
  
  ##############################################################################################################################################################################################
  #Map clusters back to genes
  ##############################################################################################################################################################################################
  
  print("Clustering using HAC")
  cluster = 5
  som_cluster <- cutree(hclust(dist(getCodes(som_model))), cluster)
  
  #Make dataframe for use later
  som.data <- as.data.frame(som_input)
  
  # get vector with cluster value for each original data sample
  cluster_assignment <- som_cluster[som_model$unit.classif]
  
  # for each of analysis, add the assignment as a column in the original data:
  sig$module <- cluster_assignment
  
  print("Returning genes with clusters")
  
  return(sig)
  
  
  
  
}



# 
# 
# 
# # print(file_name)
# # 
# # #Set which cells to be used (so as to subset large data set for comparison)
# # 
# # file_name_number <- 1
# # cell1 <- "EM"
# # cell2 <- "N"
# # cluster <- 5
# # rlen <- 2000
# # #Decide based on genes < 3000 use small > 3000 use big
# # small_big <- FALSE
# # #Select whether you want prompts etc (use when developing shiny app)
# # user_interface <- FALSE
# # 
# # 
# # genes <- som_return_cluster(file_name_number = 1, 
# #                             cell1 = "EM",
# #                             cell2 = "N", 
# #                             user_interface = FALSE, 
# #                             cluster = 5, 
# #                             rlen = 1000, 
# #                             small_big = FALSE)
# 
# 
# 
# 
# 
# 
# 
# ##############################################################################################################################################################################################
# #Visualise SOM
# ##############################################################################################################################################################################################
# 
# #Visualise the SOM model results
# #Plot of the training progress - how the node distances have stabilised over time.
# #Custom palette as per kohonen package (not compulsory)
# #Palette defined by kohonen package
# 
# coolBlueHotRed <- function(n, alpha = 1) {
#   rainbow(n, end=4/6, alpha=alpha)[n:1]
# }
# 
# clusterSom(som) 
# dev.off()
# plot(som_model, type = "changes")
# #counts within nodes
# plot(som_model, type = "counts", main="Node Counts", palette.name=coolBlueHotRed)
# #map quality
# plot(som_model, type = "quality", main="Node Quality/Distance", palette.name=coolBlueHotRed)
# #neighbour distances
# plot(som_model, type="dist.neighbours", main = "SOM neighbour distances", palette.name=grey.colors)
# #code spread
# plot(som_model, type = "codes")
# 
# #######
# #Heatmap that works
# #######
# 
# setwd("~/alec/results/Plot_genes/Central_Memory")
# for(i in 1:10){
#   jpeg(colnames(getCodes(som_model))[i]
# , width = 350, height = 350)
#   plot(som_model, type = "property", property = getCodes(som_model)[,i],
#              main=colnames(getCodes(som_model))[i], palette.name=coolBlueHotRed)
#   dev.off()
# }
# i <- 1
# par(mfrow=c(2,2))
# for(i in 1:4){
#   # Get the map of these areas and filter for Dublin areas.
#   # head(data_raw)
#   # if (small_areas){
#   #   data_raw <- read.csv("./census_data/AllThemesTablesSA.csv")  
#   #   ireland_map <- readOGR('./boundary_files/Census2011_Small_Areas_generalised20m.shp', encoding = 'utf8')
#   #   #Note that the map polygons and the census data are not in the same order - rearrangement:
#   #   data_raw <- data_raw[match(ireland_map$SMALL_AREA, data_raw$GEOGDESC),]
#   #   idcol="GEOGDESC"
#   #   
#   # } else {
#   #   data_raw <- read.csv("./census_data/AllThemesTablesED.csv")  
#   #   names(data_raw)[1] <- "GEOGID" 
#   #   ireland_map <- readOGR('./boundary_files/Census2011_Electoral_Divisions_generalised20m.shp', encoding='utf8')
#   #   ireland_map$CSOED <- paste0("E", ireland_map$CSOED)
#   #   #Note that the map polygons and the census data are not in the same order
#   #   data_raw <- data_raw[match(ireland_map$CSOED, data_raw$GEOGID),]
#   #   idcol="GEOGID"
#   # }
#   
#   # #Filter now for certain counties
#   # if (filter){
#   #   counties <- c("Fingal", "Dublin City", "South Dublin", "Dn Laoghaire-Rathdown")
#   #   plot_idx <- ireland_map$COUNTYNAME %in% counties
#   #   data_raw <- data_raw[plot_idx,]
#   #   ireland_map <- ireland_map[plot_idx,]
#   #   rm(counties, filter, plot_idx)  
#   # }
#   
#   ### -------------- Data processing -------------------------
#   
#   plot(som_model, type = "property", property = getCodes(som_model)[,i],
#        main=colnames(getCodes(som_model))[i], palette.name=coolBlueHotRed)
#   
# }
# 
# 
# dev.off()
# plot(som_model, type = "property", property = getCodes(som_model)[,i],
#      main=colnames(getCodes(som_model))[i], palette.name=coolBlueHotRed)
# 
# 
# # Plot the heatmap for a variable at scaled / normalised values
# var <- 4 #define the variable to plot
# plot(som_model, type = "property", property = getCodes(som_model)[,var], main=colnames(getCodes(som_model))[var], palette.name=coolBlueHotRed)
# 
# # Plot the original scale heatmap for a variable from the training set:
# var <- 2 #define the variable to plot
# var_unscaled <- aggregate(as.numeric(data_train[,var]), by=list(som_model$unit.classif), FUN=mean, simplify=TRUE)[,2]
# plot(som_model, type = "property", property=var_unscaled, main=names(data_train)[var], palette.name=coolBlueHotRed)
# rm(var_unscaled)
# 
# #plot a variable from the original data set (will be uncapped etc.)
# # This function produces a menu for multiple heatmaps if a factor or character is chosen
# source('plotHeatMap.R')
# 
# 
# # A menu of all variables should be displayed if variable=0 
# # (note on Mac this will required working XQuartz installation.)
# plotHeatMap(som_model, data, variable=0)
# 
# # ------------------ Clustering SOM results -------------------
# dev.off()
# # show the WCSS metric for kmeans for different clustering sizes.
# # Can be used as a "rough" indicator of the ideal number of clusters
# mydata <- getCodes(som_model)
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(mydata,
#                                      centers=i)$withinss)
# par(mar=c(5.1,4.1,4.1,2.1))
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares", main="Within cluster sum of squares (WCSS)")
# 
# # Form clusters on grid
# ## use hierarchical clustering to cluster the codebook vectors
# som_cluster <- cutree(hclust(dist(getCodes(som_model))), 3)
# length(som_cluster)
# 
# plot(som_model, type = "property", property = getCodes(som_model)[,i], bgcol = pretty_palette[som_cluster],
#      main=colnames(getCodes(som_model))[i], palette.name=coolBlueHotRed)

# # Show the map with different colours for every cluster						  
# plot(som_model, type="mapping", bgcol = pretty_palette[som_cluster], main = "Clusters")
# coolBlueHotRed(som_cluster)
# 
# add.cluster.boundaries(som_model, som_cluster)
# 
#   plot(som_model, type = "property", property = getCodes(som_model)[,i], bgcol = pretty_palette[som_cluster],
#      main=colnames(getCodes(som_model))[i], palette.name=coolBlueHotRed)
# 
# 
# #show the same plot with the codes instead of just colours
# plot(som_model, type="codes", bgcol = palpal[som_cluster], main = "Clusters")
# 
# fpalpal <- coolBlueHotRed()
# pretty_palette 
# add.cluster.boundaries(som_model, som_cluster)
# 
# 
# length(som_cluster)
# 
# palpal <- c("#0000FFFF", "#00FFEFFF", "#1FFF00FF", "#FFE500FF", "#FF0000FF")
# 
# # -------------------- MAPPING OF SMALL AREAS (GEO) --------------------------
# # Plot the map of ireland, coloured by the clusters the map to show locations.
# 
# #plotting map with ggplot requires some data preprocessing.
# 
# #get the colour for each area defined in the data
# cluster_details <- data.frame(id=data$id, cluster=som_cluster[som_model$unit.classif])
# 
# # WARNING - these operations are computationally heavy (~ 30 seconds).
# if (small_areas){
#   mappoints <- fortify(ireland_map, region="SMALL_AREA")
#   mappoints <- merge(mappoints, data, by="id")
#   mappoints <- merge(mappoints, cluster_details, by="id")  
# } else {
#   mappoints <- fortify(ireland_map, region="CSOED")
#   mappoints <- merge(mappoints, data, by="id")
#   mappoints <- merge(mappoints, cluster_details, by="id") 
# }
# rm(cluster_details)
# 
# # Finally map the areas and colour by cluster
# ggplot(mappoints) + aes(long, lat, group=group, fill=factor(cluster)) + geom_polygon()  + coord_equal() + scale_fill_manual(values = pretty_palette) + 
#   geom_path(colour="white", alpha=0.5, size=0.05) # if you want an outline
# 
# 
# 
# 
# # Get the map of these areas and filter for Dublin areas.
# # head(data_raw)
# # if (small_areas){
# #   data_raw <- read.csv("./census_data/AllThemesTablesSA.csv")  
# #   ireland_map <- readOGR('./boundary_files/Census2011_Small_Areas_generalised20m.shp', encoding = 'utf8')
# #   #Note that the map polygons and the census data are not in the same order - rearrangement:
# #   data_raw <- data_raw[match(ireland_map$SMALL_AREA, data_raw$GEOGDESC),]
# #   idcol="GEOGDESC"
# #   
# # } else {
# #   data_raw <- read.csv("./census_data/AllThemesTablesED.csv")  
# #   names(data_raw)[1] <- "GEOGID" 
# #   ireland_map <- readOGR('./boundary_files/Census2011_Electoral_Divisions_generalised20m.shp', encoding='utf8')
# #   ireland_map$CSOED <- paste0("E", ireland_map$CSOED)
# #   #Note that the map polygons and the census data are not in the same order
# #   data_raw <- data_raw[match(ireland_map$CSOED, data_raw$GEOGID),]
# #   idcol="GEOGID"
# # }
# 
# # #Filter now for certain counties
# # if (filter){
# #   counties <- c("Fingal", "Dublin City", "South Dublin", "Dn Laoghaire-Rathdown")
# #   plot_idx <- ireland_map$COUNTYNAME %in% counties
# #   data_raw <- data_raw[plot_idx,]
# #   ireland_map <- ireland_map[plot_idx,]
# #   rm(counties, filter, plot_idx)  
# # }
# 
# ### -------------- Data processing -------------------------
# 
# 
# #rlog_counts size_factor_counts variance_stabilised_counts
# 
# 
# #convert the data from summations to percentages such 
# #that the characteristics of each area will be comparable.
# # source("convertCSOdata.R")
# # data <- convertCSOdata(data_raw, idcol=idcol)
# 
# #Create SOM for Census data - simple as data is well behaved.
