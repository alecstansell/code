#Write function for distrbution plots
distributionplot <- function(counts, name){
  #Melt matrix into the right form
  x = melt(as.matrix(counts))
  #Get rid of any plots currently being visualised
  #dev.off()
  #Crete column names for melted object
  colnames(x) = c('gene_id', 'Sample', 'value')
  #Plot logged data using the number 
  ggplot(x, aes(x=value, color=Sample)) + geom_density() + ggtitle(name)
}
# #Melt matrix into the right form
# x = melt(as.matrix(raw_counts))
# #Get rid of any plots currently being visualised
# dev.off()
# #Crete column names for melted object
# colnames(x) = c('gene_id', 'sample', 'value')
# #Plot logged data using the number 
# ggplot(x, aes(x=value, color=condition)) + geom_density() + ggtitle('name')

