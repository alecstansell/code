##############################################################
##Plot sample heatmap
###############################################################
plot_sample_heatmap <- function(counts, metadata, name){
  
  #Write number of cell types or conditions
  num_conditions <- nlevels(metadata$condition)
  
  #Create object pal with colours from brewer (whichcreates the nummber of colours needed for the different conditions) color ramp then interpolates these to have a gradient in the heatmap of color
  pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
  
  #Assign colour to each respective sample
  cond_colors <- pal[as.integer(metadata$condition)]
  
  #Create a heatmap using gplots function heatmap2 # cor calculates the correlation between all the samples from the count data rowside colours for the side bar for each condition, 
  #Use margins to resize heatmap to fit without changing scaling +
  heatmap.2(cor(counts), RowSideColors=cond_colors,
            trace='none', main=name, margins=c(12,15))
  
  
}

#paste("Sample Correlations", name, sep = " ")