require(ape)
require(phytools)


#' proposal
#' Reroot the input tree by moving along the branches by step size
#' drawn from uniform (0, delta).  If this step crosses an internal
#' node, go up the left branch with probability 0.5.
#' @param phy {ape::phylo}:  input rooted tree
#' @param delta {double}:  maximum step size, s.t. proposal is U(-delta, delta)
proposal <- function(phy, delta=NA) {
  if (!is.rooted(phy)) {
    stop("Input tree must be rooted.")
  }
  if (is.na(delta)) {
    # default to 1% of entire tree length
    delta <- 0.01*sum(phy$edge.length)
  }
  if (delta > sum(phy$edge.length)) {
    # FIXME: diameter (longest path) would be better
    stop("Input delta is greater than total tree length.")
  }
  
  # process initiates at root of input tree
  root <- Ntip(phy) + 1
  step <- runif(1, 0, delta)
  
  # select a descendant branch at random
  child <- sample(phy$edge[which(phy$edge[,1]==root), 2], 1)
  clen <- phy$edge.length[which(phy$edge[,2] == child)]
  
  while (TRUE) {
    print(paste(root, child, clen, step))
    
    if (step > clen) {
      # step exceeds branch length
      root <- child
      step <- step - clen
    
      if (!is.element(child, phy$edge[,1])) {
        # child is terminal branch - reflect step for symmetric proposal
        phy <- reroot(phy, child)   # reroot at tip
        
        # update child
        prev.root <- phy$edge[1,][phy$edge[,2] == child]
        children <- phy$edge[which(phy$edge[,1]==prev.root), 2]
        child <- sample(list(children[children != child]), 1)[[1]]
      }
      else {
        # climb to next internal node, update child
        child <- sample(phy$edge[which(phy$edge[,1]==root), 2], 1)
      }
      clen <- phy$edge.length[which(phy$edge[,2] == child)]
    }
    else {
      # easiest case, step within current child branch
      return(reroot(phy, child, position=step))
    }
  }
}



#set.seed(1); phy <- rtree(5); write.tree(phy)
phy <- read.tree(text="(t2:0.06178627047,((t1:0.6870228467,(t3:0.76984142,t4:0.4976992421):0.3841037182):0.1765567525,t5:0.7176185083):0.2059745749);")

par(mfrow=c(1,2))
plot(phy); plot(proposal(phy, delta=0.5))



# TODO: calculate likelihood of rooted tree by regressing divergence on dates

# TODO: Metropolis-Hastings - from another package?
