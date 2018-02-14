source("symb_ent.R")

# Temporal ordering using symbolic transfer entropy
STEorder <- function(data, ord = 5, l = 1, distance.method = "euclidean", cluster.size) {
  # data: a matrix with cells in rows and variables in columns
  # ord: the markov order or embedding dimension
  # l:the time delay (for y -> time delay l -> x)
  # distance.method: method to calculate distance(similarity) between cells
  # cluster.size: number of cells per cluster
  
  # Main steps:
  # 1. Cluster data
  # 2. Calculate STE between clusters
  # 3. Order clusters based on STE
  
  # Some checks
  data <- as.matrix(x = data)
  
  # Calculate similary matrix
  # TODO
  # Try different similarity metrics; correlation, optimal transfer
  data.dist <- dist(x = data, method = distance.method, diag = F, upper = F)
  
  # Hierarchical clustering
  # TODO
  # Try a minimum spanning tree instead of hierarchical clustering
  hc <- hclust(d = data.dist, method = "complete")
  # Use ordering from hc to split into equal size clusters
  cells.ordered <- hc$labels[hc$order]
  cluster <- split(x = cells.ordered, f = ceiling(seq_along(cells.ordered) / cluster.size))
  # NOTE: ste not defined for unequal length vectors
  
  # Calculate STE between clusters
  # TODO
  # How to calculate ste between matrices?
  # Option 1: Calculate between individual variable vectors and then combine somehow to one value.
  # Option 2: Define how to turn matrices into symbolic order and then calculate ste.
  
  # Order clusters based on STE
  # TODO 
  # Think how to do this wisely..
  # How to define individual cell order (at this point we only have cluster order)
}