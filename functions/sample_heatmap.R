#Write function for Sample heatmap
sampleHeatmap <- function(counts, name){
  heatmap.2(cor(counts), RowSideColors=cond_colors,
            trace='none', main=name, margins=c(12,15))
  print(praise("${Exclamation}! Your heatmap is ${adjective}!"))
}
