

# Load necessary libraries
install.packages("fpc") 
install.packages("cluster")   
install.packages("factoextra") 
library(cluster)
library(factoextra)
library(fpc)


# LOAD DATA INTO R #
data <- read.table("micro_H3_vs_BC_pcs_coordinates-luke.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)



### SILHOUETTE SCORING ###
# Perform K-means clustering (for example, with 2 clusters (aka centers))
set.seed(123)  # Set seed for reproducibility
clusters <- kmeans(data, centers = 2)$cluster

# Compute the distance matrix
dist_matrix <- dist(data)

# Compute silhouette scores
silhouette_scores <- silhouette(clusters, dist_matrix)

# Print silhouette scores
print(silhouette_scores)

# Plot silhouette scores
fviz_silhouette(silhouette_scores)




### DUNN INDEX ###
set.seed(123)  # Set seed for reproducibility
clusters <- kmeans(data, centers = 2)$cluster
# Compute the distance matrix
dist_matrix <- dist(data)

# Calculate the Dunn index
dunn_index <- cluster.stats(dist_matrix, clusters)$dunn
print(dunn_index)

