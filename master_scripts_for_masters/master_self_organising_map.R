##############################################################################################################################################################################################
#Master Self Organising Map Script
##############################################################################################################################################################################################

##############################################################################################################################################################################################
#Losely based on the script by by Shane Lynn 14-01-2014
##############################################################################################################################################################################################

#Make sure results are reproducible
set.seed(1)

##############################################################################################################################################################################################
#Load Packages
##############################################################################################################################################################################################

#List packages
packages <- c("kohonen", "dummies", "ggplot2", "limma", "rgdal", "rgdal", "sp", "reshape2", "maptools", "rgeos")

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
setwd('~/data')

#Load normalised count data for use in SOM
load("final_heat_normalised.RData")

#Set specific normalisation version
counts_input <- normalised_clean$variance_stabilised_counts

#Get file names for DESeq2 differentially expressed genes
file_name <- list.files("~/pairwise_no_dups")

#Choose dataset to use
file_choice <- file_name[3]

# Import DE genes from file if you are not going to use the above
sig <- read.csv(file.path("~/pairwise_no_dups", file_choice))
sig <- sig[,-1]
colnames(sig) <- c("gene_names", "log2fold")


#Set which cells to be used (so as to subset large data set for comparison)
cell1 <- "EM"
cell2 <- "N"

#Pick cells for comparison
som_input <- counts_input[, which(grepl(cell1, colnames(counts_input)) | grepl(cell2, colnames(counts_input)) == TRUE)]


#Select whether you want prompts etc (use when developing shiny app)
user_interface <- FALSE

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
#data_train <- som_input
data_train <- counts_input
som_input <- counts_input
#Train SOM using Kohonen method
data_train_matrix <- as.matrix(scale(data_train))

?scale()

#Set names
names(data_train_matrix) <- names(data_train)

#Decide based on genes < 3000 use small > 3000 use big
small_big <- TRUE

if (small_big){
  # larger grid for the small areas example (more samples)
  som_grid <- somgrid(xdim = 25, ydim = 25, topo="hexagonal")
} else {
  som_grid <- somgrid(xdim = 10, ydim=10, topo="hexagonal")
}

# SOM Model
system.time(som_model <- som(data_train_matrix,
                             grid=som_grid,
                             rlen=5000,
                             alpha=c(0.1,0.01),
                             keep.data = TRUE ))
setwd("/home/studentsgh129/doorknobdave/alecS/alec_code/coexpression/final_coexpression/results/30157_genes")
save(som_model, file = "som_model")

##############################################################################################################################################################################################
#Map clusters back to genes
##############################################################################################################################################################################################

#Select Cluster
cluster <- 5
som_cluster <- cutree(hclust(dist(getCodes(som_model))), cluster)

#Make dataframe for use later
som.data <- as.data.frame(som_input)

# get vector with cluster value for each original data sample
cluster_assignment <- som_cluster[som_model$unit.classif]

# for each of analysis, add the assignment as a column in the original data:
som.data$module <- cluster_assignment

length(som.data$SCOPE_14952010.PBMC_N)

getmodule(som.data, 2)


sig

som.data$SCOPE_14952010.PBMC_N
som_input


##############################################################################################################################################################################################
#Visualise SOM
##############################################################################################################################################################################################

#Visualise the SOM model results
#Plot of the training progress - how the node distances have stabilised over time.
#Custom palette as per kohonen package (not compulsory)
#Palette defined by kohonen package

coolBlueHotRed <- function(n, alpha = 1) {
  rainbow(n, end=4/6, alpha=alpha)[n:1]
}

clusterSom(som)
dev.off()
plot(som_model, type = "changes")
#counts within nodes
plot(som_model, type = "counts", main="Node Counts", palette.name=coolBlueHotRed)
#map quality
plot(som_model, type = "quality", main="Node Quality/Distance", palette.name=coolBlueHotRed)
#neighbour distances
plot(som_model, type="dist.neighbours", main = "SOM neighbour distances", palette.name=grey.colors)
#code spread
plot(som_model, type = "codes")

#######
i <- 1
#Heatmap that works
#######
getwd()
setwd("/home/studentsgh129/doorknobdave/alecS/alec_code/coexpression/final_coexpression/results/30157_genes/SOM/heatmaps")
for(i in 1:24){
  jpeg(paste(colnames(getCodes(som_model))[i],".jpeg", sep = '') , width = 350, height = 350)
  plot(som_model, type = "property", property = getCodes(som_model)[,i],
             main=colnames(getCodes(som_model))[i], palette.name=coolBlueHotRed)
  dev.off()
}




plot(som_model, type = "property", property = getCodes(som_model)[,i],
     main=colnames(getCodes(som_model))[i], palette.name=coolBlueHotRed)


attach(mtcars)
par(mfrow=c(2,2))
plot(wt,mpg, main="Scatterplot of wt vs. mpg")
plot(wt,disp, main="Scatterplot of wt vs disp")
hist(wt, main="Histogram of wt")
boxplot(wt, main="Boxplot of wt")











i <- 1
p
