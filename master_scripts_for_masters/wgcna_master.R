##############################################################################################################################################################################################
#Master Coexpression Script
##############################################################################################################################################################################################

##############################################################################################################################################################################################
#Adapted and written by Alec Stansell
#Based on the tutorial written by available from https:/github.com/iscb-dc-rsg/2016-summer-workshop/tree/master/3B-Hughitt-RNASeq-Coex-Network-Analysis/tutorial
# Video lecture walking through some key concepts can be found at https:/www.youtube.com/watch?v=OdqDE5EJSlA
##############################################################################################################################################################################################

#Make sure results are reproducible
set.seed(1)


##############################################################################################################################################################################################
#Install Packages
##############################################################################################################################################################################################

#Install any packages you do not currently have

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Homo.sapiens")
# # BiocManager::install("Homo.sapiens")
#
# install.packages("Homo.sapiens")


##############################################################################################################################################################################################
#Load Packages
##############################################################################################################################################################################################

packages <- c("gplots", "tidyverse", "riskyr", "ggplot2", "knitr", "limma", "reshape2", "dplyr", "RColorBrewer", "WGCNA", "praise", "Homo.sapiens","RColorBrewer")
lapply(packages, require, character.only = TRUE)

######################################################################################################################################
#Read in data
##############################################################################################################################################################################################

metadir <- "~/"
inputdir <- "~/"
functiondir <- "~/"

}else{
  functiondir <- "~/final_coexpression"
}

#Read meta data in
samples <- read.table(file.path(metadir, "name_counts_file.csv"), header = TRUE, sep = ",", col.names = c("sample_id", "condition", "patient"), row.names = NULL, colClasses = c("character", "factor", "factor"))
#kable(samples)

#Read in count data that is summarised to gene level see doc work flow for a list of script in order
raw_counts <- read.table(file.path(inputdir, "genelevelabundance_file.csv"), sep = ",", row.names=1)


##############################################################################################################################################################################################
#Get annotation data if necessary - See first part of cluster profiler script for an alternative to getting anotation data
##############################################################################################################################################################################################

# Load Homo sapiens package to use for anotation of the gene ids
#library('Homo.sapiens')

#See different types of data available
# keytypes(Homo.sapiens)
# columns(Homo.sapiens)

#Sample query
# gene_ids <- head(keys(Homo.sapiens, keytype='ENSEMBL'), 2)
#
# fselect(Homo.sapiens, keytype='ENSEMBL', keys=gene_ids,
#        columns=c('ALIAS', 'TXCHROM', 'TXSTART', 'TXEND'))



######################################################################################################################################
#Data Prepartion
######################################################################################################################################

#Are you running this script for the first time or are you loading from saved R.dat objects?

loading_from_file <-  FALSE


#Remove potential otuliers
input <- remove_outliers_sample_tree(t(raw_counts), samples)


#Remove outliers returns a list containing adjusted counts and meta data with outliers removed
Set inputs for further normalisation
input_counts <- t(input[[1]])
input_meta <- input[[2]]



###########################################################################################################################################################################################
# Normalissation
###########################################################################################################################################################################################


#Do you want to run normalisation or load data
if(loading_from_file == FALSE){

  #Load normalisation script
  source('~/functions/normalisation_coexpression.R')

  #There are two things you can normalise - cleaned counts from the above functions and raw unprocessed counts
  #Deseq normalisation of clean counts
  normalised_clean <- normalise_it(heat_counts, heat_meta)

  #Rename
  normalised_heat <- normalised_clean

}else{
  #Set word dir
  setwd('~/data')

  #save(normalised_clean, file = "final_heat_normalised.RData")
  load("final_heat_normalised.RData")

  #Rename
  normalised_heat <- normalised_clean

}



###########################################################################################################################################################################################
# Data Clean for WGCNA post normalised
###########################################################################################################################################################################################

if(loading_from_file == FALSE){

  #Would you like to use clean counts
  print("Would you like to use clean counts?
      Type 1 for clean counts
      Type 2 to use raw unprocessed counts")
  clean <- scan(file = "", nmax = 1)
  if (is.null(clean)){
    clean <- 1

  }
  if(clean == 1){
    #Load clean function
    source(file.path(functiondir, "clean_coexpression.R"), echo = TRUE)

    #Clean counts for WGCNA processing
    clean_counts <- clean_samples(raw_counts = raw_counts)

  } else{
    clean_counts <- raw_counts
    print("You chose to use raw counts")


  }

  #Would you like to remove outliers for WGCNA
  print("Would you like to remove outliers using WGCNA?
      Type 1 for clean counts
      Type 2 to use raw unprocessed counts")
  clean <- scan(file = "", nmax = 1)

  #Remove outliers based on WGCNA
  source(file.path(functiondir, "remove_outliers.R"), echo = TRUE)
  remove_outliers_sample_tree(raw_counts, samples)


} else

if(loading_from_file == FALSE){
  #Clean counts for WGCNA processing
  clean_heat <- clean_samples(raw_counts = heat_counts)
  plot_sample_heatmap(clean_heat, heat_meta, "Cleaned Final Counts Heatmap")

}


##############################################################################################################################################################
#Visualise normalised counts
##############################################################################################################################################################

##############################################################
##Plot sample heatmap
###############################################################

#Source Heatmap Function
source('~/functions/plot_sample_heatmap.R')

names(normalised_clean) <- c("Size Factor Normalised Counts", "Variance Stabilised Normalised Counts", "Rlog Normalised Counts", "Raw Counts")
#Decide if you are saving heatmaps or not
savingnewheatmaps <-  TRUE
if(savingnewheatmaps == TRUE){

  #Plot normalised heatmaps
  setwd('~/results/30157_genes/heatmaps')

    for(i in 1:4){
    #Set up jpeg file type
    jpeg(filename = paste(names(normalised_clean)[i], "_heatmap.jpeg", sep = ""),
         width = 1000, height = 700, units = "px", pointsize = 12,
         quality = 75,
         bg = "white", res = NA,
         type = c("cairo", "Xlib", "quartz"))
    plot_sample_heatmap(normalised_clean[[i]], normalised_clean[[5]], names(normalised_clean)[i])
    #Save to null device
    dev.off()
  }

}

#Plot sample heatmap for final raw counts with outliers removed
plot_sample_heatmap(normalised_heat$variance_stabilised_counts, normalised_clean$`Input Metadata`, "Variance Stabilised Counts")

#Plot sample heatmap for inital raw counts
plot_sample_heatmap(normalised_heat$input_counts, normalised_clean$`Input Metadata`, "Raw Counts")

# Get rid of plot for the next plot
graphics.off()


######################################################################################################################################
#Distributions of the different data sets
######################################################################################################################################

#Load Distribution plot function
source('~/distribution_plot.R')

#Loop through all the samples to create the different distributions for cleaned and normalised and just normalised
#Clean Normalsied for all four comparisons

#Would you like to save distributions?
savingdistributions  <- TRUE
if(savingdistributions == TRUE){

  setwd('~/results/30157_genes/distribution_plots')
  for(i in 1:4){
    i <-  1
   #Open jpeg file to create
    jpeg(filename = paste(names(normalised_clean)[i], "_distribution_plot.jpeg", sep = ""),
         width = 900, height = 480, units = "px", pointsize = 12,
         quality = 75,
         bg = "white", res = NA,
         type = c("cairo", "Xlib", "quartz"))
    distributionplot(normalised_clean[[i]], names(normalised_clean)[i]) + xlim(1, 15)
    #Save to null device
    dev.off()
  }

    #Run for raw counts
    setwd('~/results/raw')
    #Open jpeg file to create
    jpeg(filename = "Raw Counts",
         width = 900, height = 480, units = "px", pointsize = 12,
         quality = 75,
         bg = "white", res = NA,
         type = c("cairo", "Xlib", "quartz"))
    distributionplot(raw$Counts, "Raw counts")
    #Save to null device
    dev.off()


}


#Visualise Normalised and raw counts
distributionplot(raw_counts$LIR_03.PBMC_CM ,"Variance stabilised Counts") + xlim(0, 25)

#More saving options

# setwd('~/results/raw_normalised')
#
# for(i in 1:4){
#   i = 4
#   #Open jpeg file to create
#   jpeg(filename = names(normalised_raw[i]),
#        width = 900, height = 480, units = "px", pointsize = 12,
#        quality = 75,
#        bg = "white", res = NA,
#        type = c("cairo", "Xlib", "quartz"))
#   distributionplot(normalised_raw[[i]], names(normalised_raw)[i])
#   #Save to null device
#   dev.off()
#
# }

for(i in 1:4){

  i <- 4
  #Open jpeg file to create
  jpeg(filename = paste(names(normalised_heat)[i], "_distribution_plot.jpeg", sep = ""),
       width = 900, height = 480, units = "px", pointsize = 12,
       quality = 75,
       bg = "white", res = NA,
       type = c("cairo", "Xlib", "quartz"))
  distributionplot(normalised_heat[[i]], names(normalised_heat)[i]) + xlim(0, 15)
  #Save to null device
  dev.off()
}

######################################################################################################################################
#Heatmap of Transformed / Normalised Counts
######################################################################################################################################

#Load heatmap function
source('~/sample_heatmap.R')

sampleHeatmap(normalised_heat, 'test')

#Set results directory
setwd("~/alec/results/coexpression/Sample_heatmaps_normalised")

#View Heatmaps for the 4 different normalisation methods

for (i in 1:length(normalised)){
  sampleHeatmap(normalised[[i]], names(normalised)[i])

}
#Compare to sample heatmap
sampleHeatmap(size_factor_counts, "raw_counts")

#Need to improve Png parameters to write higher res files
#Write these to file in your specified directory
outdir <- "~/doorknobdave/alecS/results/Coexpression_analysis/Gene_level_analysis"
setwd('~/results/final_normalised_heatmaps')
setwd(outdir)
for (i in 1:length(normalised_)){
  png(filename=names(normalised)[i], width = 700, height = 600, )
  sampleHeatmap(normalised[[i]], names(normalised)[i])
  dev.off()

}

###########################################################################################################################################################################################
#Decide on what type of DE analysis you would like
###########################################################################################################################################################################################
#If you would like to run differential expression analysis you can use the below script

#It is advisable to rather run a separate differential analysis experiment and start from differentially expressed files as detailed in (#Get back to raw counts from pairwise deseq file)

print("Choose what type of DE Analysis you would like to use")
invisible(readline(prompt="Press [enter] to continue"))

print("Press 1 for Deseq2, Press 2 for Limma")
typeofDEanalysis <- scan(file = "", nmax = 1)


###########################################################################################################################################################################################
#Limma Differential Expression
###########################################################################################################################################################################################
if(loading_from_file == FALSE){

  if(typeofDEanalysis = "2"){

    source('~/limma_coexpression.R')

    source('~/distribution_plot.R')

    CM_EM <- unlock_the_limma(raw_counts, samples)

    CM_N <- unlock_the_limma(raw_counts, samples)

    EM_N <- unlock_the_limma(raw_counts, samples)

    limmasig_genes <- list(CM_EM, CM_N, EM_N)

    write.csv(CM_EM, "CM_EM_coexpressed_mod_nums_genes.csv")

    setwd('~/data')
    save(limmasig_genes, file = "limmasig_genes_list.RData")

  }

  else if(loading_from_file == TRUE){

    #Choose to use normalised variance stabilised counts
    counts_input <- normalised_heat$variance_stabilised_counts
    #Load limma genes and select for input
    load("limmasig_genes.RData")

    #Select variance stabilised counts for DEGs from each comparison for WGCNA
    #For CM EM
    CM_EM <- limmacountDEGS[[1]]
    CM_EM_var_stab <- counts_input[rownames(counts_input) %in% CM_EM,]

    #For CM N
    CM_N <- limmacountDEGS[[2]]
    CM_N_var_stab  <- counts_input[rownames(counts_input) %in% CM_N,]

    #For EM N
    EM_N_var_stab  <- limmacountDEGS[[3]]
    EM_N_var_stab  <- counts_input[rownames(counts_input) %in% EM_N,]

}

}


###########################################################################################################################################################################################
#Deseq2 Differential analyis
###########################################################################################################################################################################################
typeofDEanalysis <- 1

if(typeofDEanalysis = "1"){

  #Files needed for function of deseq which does tximport and gets list of all comparisons sig genes
  samples <-  read.table(file.path(metadir, "meta_for_coexpression_30.csv"), header = TRUE, sep = ",", col.names = c("sample_id", "condition", "patient"), row.names = NULL, colClasses = c("character", "factor", "factor"))
  tx2gene <- read.csv(file.path(meta_dir, "tx2gene94.txt"))
  dir <- "~/alec/data/intermediary_data/kallisto_out"

  #Load function for deseq
  scriptdir <- "~/alec/alec_code/functions"

  #Decide whether to use a specific comparion for coexpression or all genes from all coexpression
  useallcomparions <- FALSE

  #Get list of sig DESEq genes to be used for coexpression (these are mapped back to the raw counts)
  if(useallcomparions){
    source(file.path(scriptdir, "deseq2_function.R"), echo = TRUE)
    sig <- deseqit(samples, meta_dir, tx2gene, dir)
  } else{
    source(file.path(scriptdir, "deseq_only_no_loop.R"), echo = TRUE)
    sig <- deseqonly(samples, meta_dir, tx2gene, dir)

    }

}


#####################################################################################################################################################################################################
#Get back to raw counts from pairwise deseq file
######################################################################################################################################################################################################

#Location where differentially expressed files are
file_out_deseq <- "~/degs_output"
setwd(file_out_deseq)

#If saved as an Robj
load(file.path(file_out_deseq,"robject.RData"))

##################################################################
##################################################################

##################################################################
##################################################################
'Selection on which dataset, which ids, which logfc test'


#Read in sig differntially genes
sig <- "list of genes"

cell1 <- "N"
cell2 <- "EM"

#Pick Gene Ids to be used
sig$gene_names <- sig$Row.names

#Pick log 2 fold change scores
sig$log2fold <- sig$log2FoldChange

#Choose whether you want upregulated downregulated or all differentially expressed genes
interface <- FALSE
if(interface){
  #Print choices
  print("Choose whether you would like to use downregulated, upregulated or all DE genes")
  invisible(readline(prompt="Press [enter] to continue"))
  print("Type up for upregulated, down for downregulated and all for all Differentially expressed genes")
  reg <- scan(file = "", what = "character", nmax = 1)


}else{
  reg <- 'all'

}

#If statement to select specific genes needed
if(reg == "up") {
  sig <- subset(sig, log2fold > 0)
  sig <- sig[order(sig$log2fold, decreasing = FALSE),]
  print("Upregulated genes selected")
} else if(reg == "down"){
  sig <- subset(sig, log2fold < 0)
  sig <- sig[order(sig$log2fold, decreasing = TRUE),]
  print("Downregulated genes selected")
} else if(reg == "all"){
  sig <- sig
  sig <- sig[order(sig$log2fold, decreasing = TRUE),]
  print("All genes selected")


}

#Set specific counts input data frame
counts_input <- normalised_clean$variance_stabilised_counts


#Get raw counts for the differentially expressed comparison
counts_input <- counts_input[rownames(counts_input) %in% sig$gene_names,]
counts_input <- counts_input[, which(grepl(cell1, colnames(counts_input)) | grepl(cell2, colnames(counts_input)) == TRUE)]



######################################################################################################################################
######################################################################################################################################
#START WGCNA
######################################################################################################################################
######################################################################################################################################
nrow(counts_input)
#Set wgcna input depending on comparison decided above
wgcna_input <- counts_input

#Similarity matrix function
source('~/functions/similarity_matrix.R')

#Create simularity matrix
sim_matrix <- cordist(wgcna_input)

######################################################################################################################################
#Heatmap to view initial similarity matrix - notice very large clustering - this is not informative therefore a need for adjacency matrix  where powers are raised to reduce the stringency for expressed genes
######################################################################################################################################

#make it 500
heatmap_indices <- sample(nrow(sim_matrix), 500)

heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA,
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Similarity matrix',
          density.info='none', revC=TRUE)

######################################################################################################################################
#Adjacency matrix
######################################################################################################################################

#Load Scale free independence function
source('~/functions/plot_scale_independence.R')

#Plot Scale Independence to decide powers to use
plot_scale_independence_mean_connectivity(wgcna_input)

# Construct adjacency matrix using WGCNA function - values are raised to the power 12 - this is arbitary and a calculation as to the correct power to use can be found in I_clean.R
adj_matrix <- adjacency.fromSimilarity(sim_matrix, power = 9, type='signed')


#Reading regarding power choice http:/labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html (E.g. item #5)
#https:/support.bioconductor.org/p/66101/#66195

# Delete similarity matrix to free up memory
rm(sim_matrix)

#Convert to matrix
gene_ids <- rownames(adj_matrix)
adj_matrix <- matrix(adj_matrix, nrow=nrow(adj_matrix))
rownames(adj_matrix) <- gene_ids
colnames(adj_matrix) <- gene_ids
dev.off()
dim(adj_matrix)

######################################################################################################################################
#######Heatmap for adjacency matrix now note the genes that are expreswsed at a high level ocupy a very small part of the graph you can alter the powers of adjacency matrix to affect this
######################################################################################################################################

heatmap.2(t(adj_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA,
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Adjacency matrix',
          density.info='none', revC=TRUE)


############################################################################################################################################

# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.

##########

'Topological overlap matrix alternative calculation'
#prior to hclust

TOM = TOMsimilarity(adj_matrix);
dissTOM = 1-TOM
###########


###############################
'Eigenes for hub gene analysis'
################################

#See https:/horvath.genetics.ucla.edu/html/
#CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-06-RelatingToExt.pdf

#load
library(purrr)

#Change to avoid errors
options(stringsAsFactors = FALSE)

#Because you probably have no idea what you are doing
?moduleEigengenes()

#Transpose it
transpose <- t(wgcna_input)
colnames(transpose) <- rownames(wgcna_input)

# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=30,
                                   deepSplit=TRUE)
# assign a color to each module for easier visualization and referencing
module_colors <- labels2colors(module_labels)


#Use module colours =
ncol(transpose)
length(module_colors)

#Create eigenevibe to represent each module by the first principle component eigene
datME <- moduleEigengenes(transpose,module_colors)$eigengenes
MET <- orderMEs(cbind(datME))
signif(cor(datME, use="p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
hclustdatME$

# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Module Eigengene Clustering Dendrogram")
dev.off()
#Plot of Eigen Network
tiff("eigegene.tiff",
     units="in",
     width=6,
     height=6,
     res=150)
par(mfrow=c(1, 1))
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2),
                      marHeatmap = c(3,4,1,2), cex.lab = 0.8,
                      xLabelsAngle = 90)
dev.off()



#Get kmes
datKME <- signedKME(transpose, datME, outputColumnName="")



#Filter by kmes
FilterGenes= abs(datKME$black)>.9 #abs(GS1)> .2
table(FilterGenes)

dataset <- wgcna_input
nrow(dataset)

###############################
'Hub Genes Analysis'
################################


# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=30,
                                   deepSplit=TRUE)


# assign a color to each module for easier visualization and referencing
module_colors <- labels2colors(module_labels)
dataset <- t(wgcna_input)
colnames(dataset) <- rownames(wgcna_input)

hub_genes <- chooseTopHubInEachModule(dataset, module_colors,
                                        omitColors = "grey",
                                        power = 9,
                                        type = "signed")
#Create dataframe for viewing
hub_genes <- as.data.frame(hub_genes)
print(as.data.frame(hub_genes))

#Get info about genes
hub_info <- select(Homo.sapiens, keytype='ENSEMBL', keys=hub_genes$hub_genes,
                    columns=c('GENENAME', 'SYMBOL'))

# Get rid of duplicated entries with descriptions
hub_info <- hub_info[!duplicated(hub_info$ENSEMBL),]

hub_genes <- cbind(hub_genes, hub_info)
hub_genes$ENSEMBL <- NULL
print(as.data.frame(hub_genes))

sig
counts_input
sig
getwd()

#Write to file
setwd("~/results/30157_genes/NvsEM")
write.csv(hub_genes, "Hubgenes.csv")

#Bind to retain all meta data on genes including logfcs
#hub_info <- cbind(hub_info, sig)


###############################
'Dendrogram Clustering'
################################


# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")

dissTOM

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=30,
                                   deepSplit=TRUE)

# assign a color to each module for easier visualization and referencing
module_colors <- labels2colors(module_labels)

#Plot of just cluster dendrogram
# plot(gene_tree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#      labels = FALSE, hang = 0.04);

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(module_labels)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(gene_tree, mergedColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.05,
                    addGuide = TRUE, guideHang = 0.05)


gene_info <- NULL

gene_info <- select(Homo.sapiens, keytype='ENSEMBL', keys=rownames(sig),
                    columns=c('GENENAME', 'SYMBOL'))


# Get rid of duplicated entries with descriptions
gene_info <- gene_info[!duplicated(gene_info$ENSEMBL),]

#Bind to retain all meta data on genes including logfcs
gene_info <- cbind(gene_info, sig)

#Add module colours to info file
gene_info <- cbind(gene_info, module=module_colors)



gene_info$module
# Include RGB versions of module colors for better assignment in Cytoscape
gene_info$color_rgb <- col2hex(gene_info$module)

#Laad function to create pretty heatmaps from SOM script
coolBlueHotRed <- function(n, alpha = 1) {
  rainbow(n, end=4/6, alpha=alpha)[n:1]
}

#Get order number for lgo2 values to create values for heatmap in cytoscape
sig$order = findInterval(sig$log2fold, sort(sig$log2fold))

gene_info$log2foldcolour <- coolBlueHotRed(nrow(sig))[sig$order]

gene_info$log2foldchange <- sig$log2FoldChange
sig$log2fold
gene_info$rgblog2fold <- col2hex(gene_info$log2foldcolour)

data <- gene_info
data$WGCNA_module <- module_labels
setwd("~/results/30157_genes/NvsEM")
write.csv(data, "NvsEM_coexpressed_mod_nums_genes.csv")


setwd("~/results/30157_genes/all_methods")
EMvsCM_wgcna_allmethods <- read.csv("wgcna_em_cm_final.csv")
NvsCM_wgcna_allmethods <- read.csv("NCM_wgcna_final.csv")
NvsEM_wgcna_allmethods <- read.csv("wgcna_unconnected_removed.csv")

final_wgcna_allremoved <- list(EMvsCM_wgcna_allmethods,NvsCM_wgcna_allmethods,NvsEM_wgcna_allmethods)
names(final_wgcna_allremoved) <- c("EMvsCM_wgcna_allmethods","NvsCM_wgcna_allmethods","NvsEM_wgcna_allmethods")
setwd("~/results/30157_genes/all_methods")
save(final_wgcna_allremoved, file = "final_wgcna_all_methods.RData")


######################################################################################
######################################################################################

#Plot distribution to view dataset
distributionplot(adj_matrix, 'test')

#Decide threshold value to remove values close to 0
sum(adj_matrix > 0.50)

outdir <- "~/outdirectory"
setwd(outdir)

#Load Graphml function
source('~/functions/export_network_to_graphml.R')

  g <- export_network_to_graphml(adj_matrix, filename='network_file_name.graphml',
                                   threshold=0.50, nodeAttrDataFrame=gene_info)
