#########################################################################################################
#########################################################################################################
'Differential Expression Lists Master Script'
#########################################################################################################
#########################################################################################################

'This script runs DESeq2, sleuth and limmaR differential expression analysis'

#########################################################################################################
'Preparation'
#########################################################################################################

#Load Packages
packages <- c("gplots", "ggplot2", "knitr", "limma", "reshape2", "dplyr", "RColorBrewer", "WGCNA", "praise", "Homo.sapiens","RColorBrewer")
lapply(packages, require, character.only = TRUE)

#Set directories for your experiment
out_dir <- "~/"
meta_dir <- "~/"
inputdir <- "~/count_matrix"

#Read meta data in
meta_data <- read.table(file.path(meta_dir, "name_of_meta_file.csv"), header = TRUE, sep = ",", col.names = c("sample_id", "condition", "patient"), row.names = NULL, colClasses = c("character", "factor", "factor"))

#Read in count data that is summarised to gene level see doc work flow for a list of script in order
raw_counts <- read.table(file.path(inputdir, "name_raw_count_file.csv"), sep = ",", row.names=1)

#In the case you are loading an R data file (from normalised counts etc)
setwd('')
load("normalised_counts.RData")

#Set raw and meta data
raw_counts <- #file containing columns of samples and rows of gene counts
meta_data <- #meta data describing sample name as rows and condition per sample (and any other info you may have about samples)


#########################################################################################################
'Get Raw Counts Matrix'
#########################################################################################################
#In the case you would like to get counts from tximport

require(tximport)

##Set base directories
work_dir <- "~/"
meta_dir <-  "~/"
tx2gene <- read.csv(file.path("~/meta", "transcript2genemappingfile.txt"))

#Read in sample meta data for file names
meta <- read.csv(file.path(meta_dir, "names_of_files_to_import.csv"), header = TRUE)

#Set rownmaes of comparison to sample ids
rownames(comparison) <- meta$sample_title

#Get file paths for relevant comparison
files_import <- file.path(work_dir, meta$sample, "abundance.h5")

#Name correctly
names(files_import) <- meta$sample_title

#Run tximport to get the kallisto files set txout to false to return gene level abundance
txi.object_healthy <- tximport(files_import, type = "kallisto", txOut = FALSE, tx2gene = tx2gene)

counts <- txi.object_healthy$counts
colnames(counts) <- meta$sample

#Write to file
write.csv(counts, file.path('~/','name_count_file.csv' ))

raw_counts <- counts
meta_data <- meta


#########################################################################################################
'DESeq2'
#########################################################################################################

##############################################
#Run DESeq2
##############################################

#Files needed for function of deseq which does tximport and gets list of all comparisons sig genes

#samples meta data
samples <- meta_data

#Transcript to gene mappings
tx2gene <- read.csv(file.path(meta_dir, "transcript2genemappingfile.txt"))
meta_dir <-  "~/"

dir <-"~/"

file_out_deseq <- "~/"

#Set function directory
scriptdir <- "~/functions"

#Load DE
source('~/functions/deseq_only_no_loop.R')
#source('~/functions/deseq_return_dds.R')

#Run Deseq2 returning alist of sig genes and a dds object
EMvsCM_sig_genes <- deseqonly(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
NvsCM_sig_genes <- deseqonly(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)
NvsEM_sig_genes <- deseqonly(samples = samples, meta_dir = meta_dir, tx2gene = tx2gene, dir = dir)

#Select dds object for saving purposes
ddsEMvsCM <- EMvsCM_sig_genes[[2]]
ddsNvsCM <- NvsCM_sig_genes[[2]]
ddsNvsEM <- NvsEM_sig_genes[[2]]

#Save DDS Objects
deseq_ddsobjects <- list(EMvsCM_sig_genes, NvsCM_sig_genes, NvsEM_sig_genes)
names(deseq_ddsobjects) <- c("EMvsCM_sig_genes", "NvsCM_sig_genes", "NvsEM_sig_genes")
save(deseq_ddsobjects, file = file.path(file_out_deseq,"DDSobjects_deseq.RData"))



#########################################################################################################
'LimmaR Differential Expression'
#########################################################################################################


############################################
#Load Functions
############################################

#Load Limma Script
source('~ /coexpression/final_coexpression/limma_coexpression.R')

#Load Distribution Plot Script
source('~ /coexpression/final_coexpression/distribution_plot.R')

############################################
#Run Limma using function
############################################
library(dplyr)

#Get gene lists for limmaR comparisons - in this case 3
CM_EM <- unlock_the_limma(raw_counts, meta_data)
CM_N <- unlock_the_limma(raw_counts, meta_data)
EM_N <- unlock_the_limma(raw_counts, meta_data)

#Generate List to store genes in
limmasig_latency <- list(CM_EM, CM_N, EM_N)

#Name list items for easy storage
names(limmasig_latency) <- c("CM_EM", "CM_N", "EM_N")

############################################
#Save Files
############################################

#Set dir you want files to go
file_out <- "~ /degs_output/latencystudy/limma"

#Save r object
save(limmasig_latency, file = file.path(file_out,"limmasig_genes_list.RData"))

#load(file.path(file_out,"limmasig_genes_list.RData"))

#Save all lists
for(i in 1:3){

  name <- paste(names(limmasig_latency)[i], "limma_significant_no_dups.csv", sep = "_")

  write.csv(limmasig_latency[[i]], file.path(file_out, name))


}

##############################################
#Select Upregulated Or Downregulated Genes
##############################################

sig <- deseq_gene_list$EMvsCM
sig <- deseq_gene_list$NvsCM
sig <- deseq_gene_list$NvsEM

#Choose whether you want upregulated downregulated or all differentially expressed genes
sig <- sig[,-1]
colnames(sig) <- c("gene_names", "log2fold")


#Print choices
print("Choose whether you would like to use downregulated, upregulated or all DE genes")
invisible(readline(prompt="Press [enter] to continue"))
print("Type up for upregulated, down for downregulated and all for all Differentially expressed genes")
reg <- scan(file = "", what = "character", nmax = 1)

reg <- 'up'
i <- 1


#######################################################################################################################
'Count Number of Genes upregulated and downregulated'
#######################################################################################################################

#####################################################
#Count DESeq2
#####################################################

for(i in 1:3){

  #Upregulated
  df <- deseqhealthy[[i]]
  df <- subset(df, log2FoldChange > 0)
  print(paste(names(deseqfinal8signficant[i]), nrow(df), "Upregulated genes selected", sep = " "))

  #Down
  df <- deseqhealthy[[i]]
  df <- subset(df, log2FoldChange < 0)
  print(paste(names(deseqfinal8signficant[i]), nrow(df), "Downregulated genes selected", sep = " "))



}

#####################################################
#Count Sleuth
#####################################################

nrow(list_sleuth$EMvsCM_sig_sleuth)
nrow(list_sleuth$CMvsN_sig_sleuth)
nrow(list_sleuth$EMvsN_sig_sleuth)

#####################################################
#Count Limma
#####################################################
for(i in 1:3){

  #Upregulated
  df <- limmasig_genes[[i]]
  df <- subset(df, logFC > 0)
  print(paste(names(limmasig_genes[i]), nrow(df), "Upregulated genes selected", sep = " "))

  #Down
  df <- limmasig_genes[[i]]
  df <- subset(df, logFC < 0)
  print(paste(names(limmasig_genes[i]), nrow(df), "Downregulated genes selected", sep = " "))



}

#######################################################################################################################
'Sleuth'
#######################################################################################################################

#Sleuth Code

#Use commander script for analysis
source('~/functions/commander_master_sleuth.R')



##############################################
#Venn Diagram
##############################################

#Change wd to output venn diagrams
setwd("~/degs_output/venn_diagrams")

#Load Venn Script
source('~/functions/venn_diagram.R')

#Load Deseq gene list
file_out_deseq <- "~/degs_output/HIV_control/deseq2_degs"
load(file.path(file_out_deseq,"final8significantgenes.RData"))

#Load Limma gene lists
file_out_limma <- "~/degs_output/limma_degs"
load(file.path(file_out_limma,"limmasig_genes_list.RData"))

#Load sleuth gene lists
file_out_sleuth <- "~/degs_output/sleuth_degs"
load(file.path(file_out_sleuth,"sleuthsig_genes_list.RData"))

############################################################################################
#Generate Venn Diagrams using venn script
############################################################################################

#Where you would like to produce venn diagrams
setwd('~ /degs_output/venn_diagrams')

#Venn for CM EM
vennit(list_sleuth$EMvsCM_sig_sleuth$target_id, rownames(deseqfinal8signficant$EMvsCM), limmasig_genes$CM_EM$Gene_Name, "sleuth", "DESeq2", "limma", "EMvsCM_venn.png")

#Venn for CM N
vennit(list_sleuth$CMvsN_sig_sleuth$target_id, rownames(deseqfinal8signficant$NvsCM), limmasig_genes$CM_N$Gene_Name, "sleuth", "DESeq2", "limma", "CMvsN_venn.png")

#Venn for EM N
vennit(list_sleuth$EMvsN_sig_sleuth$target_id, rownames(deseqfinal8signficant$NvsEM), limmasig_genes$EM_N$Gene_Name, "sleuth", "DESeq2", "limma", "EMvsN_venn.png")


############################################################################################
"Get Genes That Are Differentially Expressed In ALl 3 Methods"
############################################################################################

#EMvsCM All Sig For ALl Methods
final <- NULL
deseq_sleuth <- NULL
all <- NULL

#Get lists and accumulate
rownames(list_sleuth$EMvsCM_sig_sleuth) <- list_sleuth$EMvsCM_sig_sleuth$target_id
rownames(limmasig_genes$CM_EM)
deseq_sleuth <- merge(as.data.frame(deseqfinal8signficant$EMvsCM), as.data.frame(list_sleuth$EMvsCM_sig_sleuth), by="row.names", sort=FALSE)
rownames(deseq_sleuth) <- deseq_sleuth$Row.names
all <- merge(deseq_sleuth, as.data.frame(limmasig_genes$CM_EM),  by="row.names", sort=FALSE)
rownames(all) <- all$Row.names
all$Row.names <- NULL
final <- merge(all, as.data.frame(normalised_clean$variance_stabilised_counts),  by="row.names", sort=FALSE)
final$logFC <- -(final$logFC)
head(final)
EMvsCM_all <- final


#CMvsN All Sig For ALl Methods
final <- NULL
deseq_sleuth <- NULL
all <- NULL

#Get lists and accumulate
rownames(list_sleuth$CMvsN_sig_sleuth) <- list_sleuth$CMvsN_sig_sleuth$target_id
rownames(limmasig_genes$CM_N)
deseq_sleuth <- merge(as.data.frame(deseqfinal8signficant$NvsCM), as.data.frame(list_sleuth$CMvsN_sig_sleuth), by="row.names", sort=FALSE)
rownames(deseq_sleuth) <- deseq_sleuth$Row.names
all <- merge(deseq_sleuth, as.data.frame(limmasig_genes$CM_N),  by="row.names", sort=FALSE)
rownames(all) <- all$Row.names
all$Row.names <- NULL
final <- merge(all, as.data.frame(normalised_clean$variance_stabilised_counts),  by="row.names", sort=FALSE)
final$logFC <- -(final$logFC)
NvsCM_all <- final


#EMvsCM All Sig For ALl Methods
final <- NULL
deseq_sleuth <- NULL
all <- NULL

#Get lists and accumulate
rownames(list_sleuth$EMvsN_sig_sleuth) <- list_sleuth$EMvsN_sig_sleuth$target_id

deseq_sleuth <- merge(as.data.frame(deseqfinal8signficant$NvsEM), as.data.frame(list_sleuth$EMvsN_sig_sleuth), by="row.names", sort=FALSE)
rownames(deseq_sleuth) <- deseq_sleuth$Row.names
all <- merge(deseq_sleuth, as.data.frame(limmasig_genes$EM_N),  by="row.names", sort=FALSE)
rownames(all) <- all$Row.names
all$Row.names <- NULL
final <- merge(all, as.data.frame(normalised_clean$variance_stabilised_counts),  by="row.names", sort=FALSE)
final$logFC <- -(final$logFC)
head(final)
NvsEM_all <- final

#################Save results#####################

results_all_methods <- list(EMvsCM_all, NvsCM_all, NvsEM_all)
names(results_all_methods) <- c("EMvsCM_all", "NvsCM_all", "NvsEM_all")
setwd("~ /degs_output/occuring_all")
save(results_all_methods, file = "list_genes_occuring_in_all_methods.RData")
load("list_genes_occuring_in_all_methods.RData")

file_out <- "~ /degs_output/occuring_all"
#Save all lists
for(i in 1:3){

  name <- paste(names(results_all_methods)[i], "all_3_methods.csv", sep = "_")

  write.csv(results_all_methods[[i]], file.path(file_out, name))


}
