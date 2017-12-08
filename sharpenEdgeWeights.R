

##############################################################################################
# implementation of "edge sharpening" algorithm from "On Clustering Using Random Walks"      #
# by Harel and Koren (http://www.wisdom.weizmann.ac.il/~harel/papers/Clustering_FSTTCS.pdf)  #
#                                                                                            #
# version 1.0                                                                                #
##############################################################################################

require(network)
require(matrixcalc)

#### inputs ####
# G: object of class 'network'
# k: maximum walk length
# steps: steps for algorithm
# plot: plot weight histogram

sharpenEdgeWeights <- function(G, k = 3, steps = 4, plot = TRUE)
{
  #### internal functions #####
  # vector of transition probabilities
  # ith entry is transition probability from node v to node i
  ###      p(i,j) = w(i,j)/d(i)     ###
  
  prob_vector <- function(G, v)
  {
    p <- rep(0, get.network.attribute(G, "n"))
    nbhd <- get.neighborhood(G,v)
    eid <- unlist(get.dyads.eids(G, rep(v, length(nbhd)), nbhd))
    w <- get.edge.attribute(G, "weight")[eid]
    p[nbhd] <- w / sum(w)
    return(p)
  }
  
  #### similarity function ####
  ## inputs:
  # k: maxiumn walk length
  # Pk: visit matrix
  # r: row of edge list fed from apply()
  similarity1 <- function(r, k, Pk)
  {
    L1norm <- function(x) sum(abs(x))
    return(exp(2*k - L1norm(Pk[,r[1]] - Pk[,r[2]])) - 1)
  }
  

  #### main loop ####
  EL<-as.matrix.network.edgelist(G)
  n<-get.network.attribute(G,"n")

  for (i in 1:steps)
  {
    # create transition matrix
    M <- sapply(X = 1:n, FUN = prob_vector, G = G)
    # create visit matrix
    Pk <- sapply(1:k, FUN = matrix.power, x = M)
    dim(Pk) <- c(n,n,k)
    Pk <- apply(X = Pk, FUN = sum, MARGIN = c(1,2))
    # compute edge weights
    W <- apply(X = EL, FUN = similarity1, MARGIN = 1, k = k, Pk = Pk)
    set.edge.attribute(G, "weight", W)
  }
  ##################
  
  #### plot weight histogram ####
  if (plot == TRUE)
  {
    hist(W, breaks = 20, main = "Histogram of edge weights", xlab = "edge weights", probability = TRUE)
  }
  
  return(structure(list(graph = G, walkMatrix = Pk , walk_length = k, steps = steps, weight_range = range(W))))
}


