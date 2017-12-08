
# example data from "On Clustering Using Random Walks" by Harel and Koren

library(network)
library(matrixcalc)

######################
####  Example 1   ####
######################

# load data 
test_graph1 <- read.csv("./examples/test_graph1.csv", stringsAsFactors = FALSE)

# create weighted graph
G1 <- network(test_graph1, matrix.type = "edgelist", directed = FALSE)
set.edge.attribute(G1, "weight", test_graph1$weight)

# view graph
plot(G1)

# apply edge sharpening algorithm
results <- sharpenEdgeWeights(G1, k = 3, steps = 4, plot = TRUE)

# find clusters using edge deletion method
new_G1 <- findClusters(results, method = "deletion", cutoff = 0.97, plot = TRUE)

# view graph clusters
plot(new_G1$graph)


######################
####  Example 2   ####
######################

# load data 
test_graph2 <- read.csv("./examples/test_graph2.csv", stringsAsFactors = FALSE)

# create weighted graph
G2 <- network(test_graph2, matrix.type = "edgelist", directed = FALSE)
set.edge.attribute(G2, "weight", test_graph2$weight)

# view graph
plot(G2)

# apply edge sharpening algorithm
results <- sharpenEdgeWeights(G2, plot = FALSE)

# find clusters using single linkage
results_1 <- findClusters(results, method = "single", plot = TRUE)

# find clusters using normalized edge sum
results_2 <- findClusters(results, method = "normalized", plot = TRUE)
