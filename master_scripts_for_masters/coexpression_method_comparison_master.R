##############################################################################################################################################################################################
#Master Clustering Comparison Experiment Using Cluster Profiler GO Enrichment
##############################################################################################################################################################################################

##############################################################################################################################################################################################
#Load Packages and script dirs
##############################################################################################################################################################################################

#List packages
packages <- c("tidyverse", "clusterProfiler", "org.Hs.eg.db", "biomaRt", "tidyr", "forcats", "ggplot2", "praise"
              )

#Load packages
lapply(packages, require, character.only = TRUE)

#Directory where function scripts are
scriptdir <- "~/functions"


functions <- c('getmodule.R','get_entrez_function.R','clusterit.R', 'mixedtofloat.R')

for (i in 1:length(functions)){

  source(file.path(scriptdir, functions[i]), echo = TRUE)


}

#Load SOM clustering function (Returns sig list of genes with comparison)
source('~/clustering/selforganisingmap/subset_self_organising_map.R')


#Load if you want to start from adjusted data
setwd("~/coexpression/final_coexpression/data")
load("final_heat_normalised.RData")

counts_input <- normalised_clean$variance_stabilised_counts



#####################################################################################################################################
#Run Self Organising map and return coexpression module groups
#####################################################################################################################################

#Run Self Organising Map for each comparison to return
EMvCM_SOM <-  som_return_cluster(user_interface = FALSE, cluster = 5 , rlen = 2000)
NvsCM_SOM <-  som_return_cluster(user_interface = FALSE, cluster = 5, rlen = 2000)
NvsEM_SOM <-  som_return_cluster(user_interface = FALSE, cluster = 5, rlen = 2000)

#Create list of above results
SOMresults <- list(EMvCM_SOM, NvsCM_SOM, NvsEM_SOM)

setwd('~/clustering/selforganisingmap/data')
list.files()
load("SOM_pairwise_results.RData")


#####################################################################################################################################
#Kmeans analysis
#####################################################################################################################################

#Returns plot of kmeans and alist with the modules
source('~/clustering/kmeans/kmeans_all.R')

setwd('~/clustering/kmeans/data')

load_from_file = TRUE

if(load_from_file){

  #load R data
  load("kmeans_input.RData")

  #Break up list
  EMvCM_kmeans <- kmeans_mean_input$EMvsCM
  NvsCM_kmeans <- kmeans_mean_input$NvsCM
  NvsEM_kmeans <- kmeans_mean_input$NvsEM

}else{

  #Load kmeans input function which takes all relevant genes for a compariosn and averages them
  source('~/clustering/kmeans/kmeans_input.R')

  #Create average foiles for kmeans input
  EMvCM_kmeans <-  kmeans_input(counts_input)
  NvsCM_kmeans <- kmeans_input(counts_input)
  NvsEM_kmeans <- kmeans_input(counts_input)

  #Save intermediary file
  kmeans_mean_input <- list(EMvCM_kmeans, NvsCM_kmeans, NvsEM_kmeans)
  names(kmeans_mean_input) <- c("EMvsCM", "NvsCM", "NvsEM")
  save(kmeans_mean_input, file = "kmeans_input.RData")


}
nrow(NvsEM_kmeans)
load_from_file <- TRUE
if(load_from_file){

  setwd('~/clustering/kmeans/data')
  load("kmeans_module_results.RData")
  EMvCMkmeans_mod <- kmeans_mod_res$EMvCMkmeans_mod
  NvsCM_kmeans_mod <- kmeans_mod_res$NvsCM_kmeans_mod
  NvsEM_kmeans_mod <- kmeans_mod_res$NvsEM_kmeans_mod

}else{

  #RUn kmeans clustering and return list of genes
  EMvCMkmeans_mod <- kmeansit(EMvCM_kmeans, "EMvCM_kmeans")
  length(NvsCM_kmeans_mod$log2FoldChange)
  NvsCM_kmeans_mod <- kmeansit(NvsCM_kmeans, "NvsCM_kmeans")
  #NvsEM_kmeans_mod <- kmeansit(NvsEM_kmeans, "NvsEM_kmeans")
  NvsEM_kmeans_mod <- kmeansit(NvsEM_kmeans, "NvsEM_kmeans")


  setwd('~/clustering/kmeans/data')
  kmeans_mod_res <- list(EMvCMkmeans_mod, NvsCM_kmeans_mod, NvsEM_kmeans_mod)
  names(kmeans_mod_res) <- c("EMvCMkmeans_mod", "NvsCM_kmeans_mod", "NvsEM_kmeans_mod")
  save(kmeans_mod_res, file = "kmeans_module_results.RData")



}


#####################################################################################################################################
#Loop through all SOM Results and get cluster profiler gene modules
#####################################################################################################################################
#na.ommit(unique(gene_list$module)[1:3])
names(kmeans_mod_res)

SOMresults <- kmeans_mod_res


setwd("~/coexpression/final_coexpression/results/30157_genes/wgcna_lists")

EMCMwgcna <- read.csv("EMvsCM_coexpressed_mod_nums_genes.csv", stringsAsFactors = TRUE)
EMCMwgcna$modulenumber <- EMCMwgcna$WGCNA_module

NCMwgcna <- read.csv("NvsCM_coexpressed_mod_nums_genes.csv", stringsAsFactors = TRUE)
NCMwgcna$modulenumber <- NCMwgcna$WGCNA_module

NEMwgcna <- read.csv("NvsEM_coexpressed_mod_nums_genes.csv", stringsAsFactors = TRUE)
NEMwgcna$modulenumber <- NEMwgcna$WGCNA_module

head(NEMwgcna)
wgcna_all <- list(EMCMwgcna, NCMwgcna, NEMwgcna)
names(wgcna_all) <- c("EMCMwgcna", "NCMwgcna", "NEMwgcna")
save(wgcna_all, file = "wgcna_listsall.RData")
load("wgcna_listsall.RData")
gene_list <- wgcna_all$EMCMwgcna

##################################################################################################################################
##################################################################################################################################

'Loop for cluster Profiler using an object produced with all 3 comparisons WGCNA'

##################################################################################################################################
#####################################################################################################################################
#Loop for clusterprof for a particular comp
list_GO <- NULL
i <-
library(Homo.sapiens)
for(i in 1:3){

  #Set gene list input
  gene_list <- final_wgcna_allremoved[[i]]

  #Remove na values
  range <- na.omit(unique(gene_list$module))

  #Ensure correct class type for lists (ie not intergers - must be characters)
  outdir <- paste("~/coexpression/final_coexpression/results/30157_genes/cluster_profiler/all_methods/", names(final_wgcna_allremoved[i]), sep = "")
  setwd(outdir)

  for(module in range[1:length(range)]){
    genes_module <- getmodule(gene_list, module)
    #genes_entrez <- getentrez(genes_module$ENSEMBL)

    genes_entrez <- select(Homo.sapiens, keytype='ENSEMBL', keys=as.character(genes),
                        columns=c('GENENAME', 'SYMBOL',"ENTREZID"))

    # Get rid of duplicated entries with descriptions
    genes_entrez <- gene_info[!duplicated(gene_info$ENSEMBL),]

        #Inputs required list of entrez genes, universe and type of ontology required (BP, CC, MF, ALL) then writes to file each module GO terms in current working directory.
    GO_termobj <- clusterit(genes_GO = genes_entrez$ENTREZID, universe_entrez = universe_entrez, ont = 'all')
    list_GO <- list(list_GO, GO_termobj )
    print(praise("${Exclamation}! ${adjective}!"))

  }



}

##################################################################################################################################
##################################################################################################################################

'Using just a list of genes inputted from WGCNA modules'

##################################################################################################################################
##################################################################################################################################

#Loop for clusterprof for a particular comp
list_GO <- NULL
library(Homo.sapiens)

#Set gene list input
gene_list <- EMvsN_deseq2

#Remove na values
range <- na.omit(unique(gene_list$module))

#Ensure correct class type for lists (ie not intergers - must be characters)
# outdir <- paste("~/coexpression/final_coexpression/results/30157_genes/NvsEM/", names(final_wgcna_allremoved[3]), sep = "")
# setwd(outdir)
# getwd()
for(module in range[1:length(range)]){

  genes_module <- getmodule(gene_list, module)
  #genes_entrez <- getentrez(genes_module$ENSEMBL)

  genes_entrez <- select(Homo.sapiens, keytype='ENSEMBL', keys=as.character(genes_module$ENSEMBL),
                         columns=c('GENENAME', 'SYMBOL',"ENTREZID"))

  # Get rid of duplicated entries with descriptions
  genes_entrez <- genes_entrez[!duplicated(genes_entrez$ENSEMBL),]



  #Inputs required list of entrez genes, universe and type of ontology required (BP, CC, MF, ALL) then writes to file each module GO terms in current working directory.
  GO_termobj <- clusterit(genes_GO = genes_entrez$ENTREZID, universe_entrez = universe_entrez, ont = 'all')
  #list_GO <- list(list_GO, GO_termobj )
  print(paste(praise("${Exclamation}! ${adjective}!"), module, "is done"))

}


##################################################################################################################################
##################################################################################################################################

'Amalgamate all GO terms'

##################################################################################################################################
#####################################################################################################################################

source('~/functions/version_nums_data_type_conversions/almagamte_go_terms.R')

#Run through all comparisons

#name <- "EMvsCM_wgcna_allmethods"

#name <-"NvsCM_wgcna_allmethods"

name <-"NvsEM_wgcna_allmethods"

name <-"module"


file_location <- paste("~/coexpression/final_coexpression/results/30157_genes/NvsEM/cluster_profiler/", name, sep = "")

out_dir <- "~/coexpression/final_coexpression/results/30157_genes/NvsEM"

almalgamate_GO(file_location, outdir)



#####################################################################################################################################
#Bar Graph to compare GO terms

#####################################################################################################################################

require(ggpubr)

#Containing all three methods for comparison with all GO terms across all three comparisons
all <- rbind(km_all, SOM_all, wgcna_all)

#Make gene ratios percent values
all$GeneRatio <- all$GeneRatio*100

# Compare methods
ggboxplot(all, x = "Type", y = "GeneRatio", color = "type",
          add = "jitter", legend = "none") +
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(all$GeneRatio), linetype = 2)+ # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 60)+        # Add global annova p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.")                      # Pairwise comparison against all

#####################################################################################################################################
#Code to see number of modules and frequency of each module
#####################################################################################################################################

#View table of frequencies of the modules
as.data.frame(table(gene_list$module))

##Write to file
write.csv(as.data.frame(table(gene_list$module)), "frequency_deseq_modules.csv")
