
#################################
#Cluster Heatmap
#################################

som_cluster <- cutree(hclust(dist(getCodes(som_model))), 5)
setwd("~/doorknobdave/alecS/alec_code/selforganisingmap/results/EMvN")

for(i in 1:18){
  som_cluster <- cutree(hclust(dist(getCodes(som_model))), 3)
  jpeg(paste((colnames(getCodes(som_model))[i]), ".jpeg", sep = "")
       , width = 350, height = 350)
  plot(som_model, type = "property", property = getCodes(som_model)[,i], bgcol = pretty_palette[som_cluster],
       main=colnames(getCodes(som_model))[i], palette.name=coolBlueHotRed)
  add.cluster.boundaries(som_model, som_cluster)
  dev.off()
  
}
i = 1
#################################
################################3

som_model
