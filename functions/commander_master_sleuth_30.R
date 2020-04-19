
###############################################################
###############################################################


##### Sleuth Differential Analysis Pipeline #####



###############################################################
###############################################################
#Contents

# 1.    Install
# 2.    Preparation
# 3.    Sleuth experiment


###############################################################
#Very Useful information used to generate this script
###############################################################

'Tutorials: - https://hbctraining.github.io/DGE_workshop_salmon/lessons/09_sleuth.html'

'Blog post describing sleuth - https://liorpachter.wordpress.com/2015/08/17/a-sleuth-for-rna-seq/'

'DESeq vs sleuth more or less https://www.r-bloggers.com/count-data-to-log-or-not-to-log/'

###############################################################
#Install Packages
###############################################################


##Install Packages

# install.packages("tidyverse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")
# install.packages("devtools", repos = "http://cran.us.r-project.org")
# set_config(config(ssl_verifypeer = 0L))
# devtools::install_github("pachterlab/sleuth")
# install.packages("cowplot")
# install.packages("rlang")
# install.packages("Rcpp")
# install.packages("tidyverse")
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# install.packages("ggplot2")
# install.packages(‘usethis’, ‘ggplot2’, ‘purrr’, ‘tidyselect’, ‘dplyr’, ‘tibble’, ‘pkgload’, ‘pillar’)
# install.packages("purr")

###############################################################
#Set up
###############################################################

#Plot density plot of date
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
                                                     "sample"), offset = 1)

dev.off()

##Load Required Packages
library("sleuth")
library("httr")
library("cowplot")
library("rhdf5")
library("rlang")

#Set optimum number of cores for experiment
options(mc.cores = 4L)

##Set base directories
work_dir <- "~/Healthy_Cells/abundance"
meta_dir <-  "~/Healthy_Cells/meta"

setwd(work_dir)

##Set Working Directory to where all meta data files and kallisto output data is
setwd("~/doorknobdave/final/intermediary_data/kallisto_out")

##Set base directories
work_dir <- "~/doorknobdave/final/intermediary_data/kallisto_out"
meta_dir <- "~/doorknobdave/alecS/meta"

###############################################################
#Sample to covariate construction
###############################################################

#Set inputs
raw_counts <-
meta_data <-


###############################################################
#Generate an summary meta data file including info relevant only for each comparison of interest
###############################################################
meta_data$cell <- meta$condition
meta_data$sample <- meta_data$sample_id


#EMvCM
EMvCM <- meta_data[ which(meta_data$cell=='EM'|meta_data$cell=='CM'), ]
EMvCM$run_accession <- NULL
rownames(EMvCM) <- EMvCM$sample

#CMvsN
CMvsN <- meta_data[ which(meta_data$cell=='CM'|meta_data$cell=='N'), ]
CMvsN$run_accession <- NULL
rownames(CMvsN) <- CMvsN$sample

#EMvN
EMvN <- meta_data[ which(meta_data$cell=='EM'|meta_data$cell=='N'), ]
EMvN$run_accession <- NULL
rownames(EMvN) <- EMvN$sample


##############################################################
#Select relevant info for experiement
#############################################################

s2c <- meta_data

s2c$sample_id <- meta_data$sample

#Read all meta data in table consisting of everything known about samples as an object s2c - sample to covariate matrix

#EMvCM
#Add column for the individual - by removing extra info about the sample id  - this can be used later for normalisation
s2c$individual <- gsub("-.*","", s2c$sample)
s2c$individual <- as.factor(s2c$individual)
s2c$cell <- factor(s2c$cell, c('N','EM','CM'), ordered = FALSE)


#Specify file paths for each of the samples h5 abundance estimate files using the sample title names (the same as the output file names created in kallisto)
#files <- file.path(work_dir, s2c$sample, "abundance.h5")
files <- file.path(work_dir, s2c$sample_id, "abundance.h5")

#Check that you have all the abundance.h5 files produced in kallisto for all samples in the original experiment

all(file.exists(files))

#######################################

##List of all sample ids
s2c$sample <- s2c$sample_id

#convert from a factor to a character and then to a numeric for the pca plots
s2c$Plasma_HIV_RNA_.copies.mLd. <-  as.numeric(as.character(s2c$Plasma_HIV_RNA_.copies.mLd.))

##Add paths to s2c dataframe set as paths!!
s2c <- dplyr::mutate(s2c, path = files)

# s2c

##Add gene names using biomart

# mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
#                          dataset = "hsapiens_gene_ensembl",
#                          host = 'ensembl.org')
# t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
#                                      "external_gene_name"), mart = mart)
# t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
#                      ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
#
#
# # Add to sleuth object
#
# so <- sleuth_prep(s2c, target_mapping = t2g)

#Load sleuth
library(sleuth)

#Create sleuth object with object cell for comparison of cell types
need_ttg_file <- FALSE
if(need_ttg_file) {

  #Load function to generate transcript to gene mapping
  source('~/doorknobdave/alecS/alec_code/functions/get_gene_file_from_transcripts_biomart.R')

  #Run function
  t2g <- gett2g()

}

library(sleuth)
#Initialise so
so <- NULL

s2c$sample <- s2c$sample_title
t2g <- read.csv(file.path("~/doorknobdave/alecS/meta", "tx2gene94.txt"))
t2g$target_id <- t2g$transcript
s2c$sample <- s2c$sample_id


#Prepare sleuth object with relevant parameters for experiment
so <- sleuth_prep(s2c,
                  full_model = ~ condition, #model design - this case cell
                  target_mapping = t2g, #t2g file with relevant transcript gene mappings
                  read_bootstrap_tpm = TRUE, #Use bootstraps
                  aggregation_column = 'ens_gene', #Agregate gene results - to give results
                  extra_bootstrap_summary = TRUE,
                  transformation_function = function(x) log2(x + 0.5))

# so <- sleuth_prep(s2c, ~ cell, target_mapping = t2g,
#                   aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)


#Fit reduced module with factors you do not want to affect the experiment
so <- sleuth_fit(so, ~patient, 'reduced')

so <- sleuth_fit(so, ~cell + patient, 'full')

# so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_lrt(so, 'reduced', 'full')

#Run sleuth live
sleuth_live(so)
