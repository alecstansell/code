######################################################################################################################################
#Create simularity matrix
######################################################################################################################################

#'Similarity measure which combines elements from Pearson correlation and Euclidean distance.
#Funtion cordistant to create a distance matrix

cordist <- function(dat) {
  cor_matrix  <- cor(t(dat))
  
  dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
  
  sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}

?heatmap.2
