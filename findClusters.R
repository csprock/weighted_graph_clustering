
##### clustering options #####
# edge deletion: deletes with weights below a cutoff threshold. Clusters are connected components of resulting graph.
# single link: merge two adjacent nodes with highest edge weight.
# normalized sum: the distance between two clusters measured by the normalized sum of edge weights between clusters.

#### input #####: 
# es_object: the output of the sharpenEdgeWeights() function

findClusters <- function(es_object, method = c("deletion", "single","normalized"), cutoff = NULL, plot = FALSE)
{
  
  ##################### clustering functions ##########################
  
  single_link_HAC <- function(G)
  {
    
    d <- as.matrix(G, attrname = "weight")
    dMAX <- max(d)
    
    N <- network.size(G)
    
    clusters <- -(1:N)
    attributes(clusters) <- list(names = get.vertex.attribute(G, "vertex.names"))
    merge_output <- matrix(0, nrow = N-1, ncol = 2)
    heights <- vector()
    
    for (j in 1:(N-1))
    {
      # set height to largest normalized distance
      dmax <- max(d)
      heights[j] <- (dMAX - dmax)/dMAX 
      
      # get row and column indices of min
      i <- which(d == dmax, arr.ind = TRUE)
      i <- i[1, , drop = TRUE]
      
      # following R convention (from bwlewis)
      merge_output[j,] <- clusters[i][ order(clusters[i]) ]
      
      # update clusters (merge min row/column with rest of vertices already in cluster (as given by positive labels))
      clusters[ c(i, which(clusters %in% clusters[ i[clusters[i] > 0] ])) ] <- j
      
      r <- apply(d[i,], 2, max)
      d[min(i),] = d[,min(i)] = r
      d[min(i), min(i)] <- 0
      d[max(i),] = d[,max(i)] = 0
    }
    
    temp <- structure(list(merge = merge_output, height = heights, order = iorder(merge_output),
                           labels = names(clusters), method = "single", dist.method = NULL), class = "hclust")
    
    return(temp)
  }
  
  normalized_sum_HAC <- function(G)
  {
    # extract weight matrix from graph
    W <- as.matrix(G, attrname = "weight")
    
    cluster_distance <- function(r, W, clusters)
    {
      cluster_indices <- list(C1 = which(clusters == r[1]), C2 = which(clusters == r[2]))
      return( normalized_weight(W, cluster_indices))
    }
    
    normalized_weight <- function(W, cluster_indices)
    {
      if (length(intersect(cluster_indices$C1,cluster_indices$C2)) > 0) stop("Warning! Non-empty cluster intersection!")
      num <- sum(W[cluster_indices$C1, cluster_indices$C2])
      denom <- length(cluster_indices$C1)^0.5 + length(cluster_indices$C2)^0.5
      return(num/denom)
    }
    
    # initializes distance matrix
    cluster_distance_matrix <- function(W, clusters)
    {
      C <- unique(clusters)
      is <- combn(C,2)
      ds <- combn(C, 2, FUN = cluster_distance, W = W, clusters = clusters, simplify = TRUE)
      return(list(clusters = is, distances = ds))
    }
    
    
    N <- network.size(G)
    clusters <- -(1:N)
    attributes(clusters) <- list(names = get.vertex.attribute(G, "vertex.names"))
    d <- cluster_distance_matrix(W, clusters)
    
    dists <- d$distances
    pairlist <- d$clusters
    
    merge_output <- matrix(0, nrow = N - 1, ncol = 2)
    heights <- vector()
    
    for (j in 1:(N-1))
    {
      i <- which.max(dists)
      p <- pairlist[,i]
      
      # merge output (R convention)
      merge_output[j,] <- p[order(p)]
      
      # compute dendrogram height
      heights[j] <- sum(clusters == p[1])*sum(clusters == p[2])
      # update clusters
      clusters[ which(clusters == p[1] | clusters == p[2]) ] <- j
      
      # compute indices containing merged clusters for deletion
      del <- union(which(pairlist[1,] %in% p), which(pairlist[2,] %in% p))
      
      # delete merged cluster distances and pairs
      pairlist <- pairlist[,-del]
      dists <- dists[-del]
      
      # compute new distances to merged cluster
      new_pairs <- rbind(setdiff(unique(clusters),j), rep(j, N - j - 1))
      new_dists <- apply(X = new_pairs, MARGIN = 2, FUN = cluster_distance, W = W, clusters = clusters)
      
      # update distances and cluster pairs
      dists <- append(dists, new_dists)
      pairlist <- cbind(pairlist, new_pairs)
    }
    
    temp <- structure(list(merge = merge_output, height = sort(heights), order = iorder(merge_output),
                           labels = names(clusters), method = "single", dist.method = "normalized edge sum"), class = "hclust")
    return(temp)
  }
  
  edge_deletion_clustering <- function(G, cutoff)
  {
    
    if (missing(cutoff)) stop("Must specify cutoff value.")
    
    # find connected components using DFS
    find_components <- function(G)
    {
      DFSutil <- function(v)
      {
        if (visit_table$visited[v] == 0)
        {
          visit_table$visited[v] <<- 1
          #print(i)
          components<<-append(components,i)
          for (u in get.neighborhood(G,v))
          {
            DFSutil(u)
          }
        } 
      }
      
      n <- get.network.attribute(G,"n")
      visit_table <- data.frame(cbind(1:n, rep(0,n)))
      names(visit_table)<-c("vertex", "visited")
      components <- vector()
      i <- 0
      
      for (v in 1:n)
      {
        if (visit_table$visited[v] == 0)
        {
          DFSutil(v)
          i<-i+1
        }
      }
      return(components)
    }
    
    delete.edges(G, which(get.edge.attribute(G, "weight") < cutoff))
    components <- find_components(G)
    set.vertex.attribute(G, "component", components)
    
    out <- list(graph = G, clusters = components)
    return(out)
  }
  
  # orders output of merge path for plotting (courtesy of bwlewis https://github.com/bwlewis/hclust_in_R)
  iorder = function(m)
  {
    N = nrow(m) + 1
    iorder = rep(0,N)
    iorder[1] = m[N-1,1]
    iorder[2] = m[N-1,2]
    loc = 2
    for(i in seq(N-2,1))
    {
      for(j in seq(1,loc))
      {
        if(iorder[j] == i)
        {
          iorder[j] = m[i,1]
          if(j==loc)
          {
            loc = loc + 1
            iorder[loc] = m[i,2]
          } else
          {
            loc = loc + 1
            for(k in seq(loc, j+2)) iorder[k] = iorder[k-1]
            iorder[j+1] = m[i,2]
          }
        }
      }
    }
    -iorder
  }
  
  #################### input logic ###################################
  
  if (method == "deletion")
  {

    temp <- edge_deletion_clustering(es_object$graph, cutoff)
    if (plot == TRUE) plot(temp$graph)
    return(temp)
    
  } else if (method == "single") {
    
    temp <- single_link_HAC(es_object$graph)
    if (plot == TRUE) plot(temp, main = "Dendrogram, Single-Linkage")
    return(temp)
    
  } else if (method == "normalized") {
    
    temp <- normalized_sum_HAC(es_object$graph)
    if (plot == TRUE) plot(temp, main = "Dendrogram, Normalized Edge Sum")
    return(temp)

  }
}



