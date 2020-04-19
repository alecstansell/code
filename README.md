# code repository

coexpression, differential and cluster analysis used in my masters

### contents

'masters_scripts_for_masters'
  scripts that direct the core parts of my project
  they rely on functions and scripts within 'functions'

'functions'
  all functions needed

#### masterscripts

coexpression_method_comparison_master.R
  Self-organising map (with heirachical clustering), k-means and WGCNA comparison and visualisation

wgcna_master
  wgcna analysis and visualisation
  includes: hub analysis, sample heatmaps, distribution plots, outlier removal, normalisation and data cleaning

differential_expresion_master.R
  differential expression analysis with three different methods (DESeq2, limmaR and sleuth)

selforganisingmap_master.R
  Code for the full transcriptional profiling and heatmap generation for SOM
