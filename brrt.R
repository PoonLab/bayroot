require(ape)
require(phytools)

#set.seed(1); phy <- rtree(5); write.tree(phy)
phy <- read.tree(text="(t2:0.06178627047,((t1:0.6870228467,(t3:0.76984142,t4:0.4976992421):0.3841037182):0.1765567525,t5:0.7176185083):0.2059745749);")

par(mfrow=c(1,2))
plot(phy)


#' proposal
#' Reroot the input tree by moving along the branches by step size
#' drawn from uniform (0, delta).  If this step crosses an internal
#' node, go up the left branch with probability 0.5.
#' @param phy {ape::phylo}:  input rooted tree
#' @param delta {double}:  maximum step size, s.t. proposal is U(-delta, delta)
proposal <- function(phy, delta=0.1) {
  if (!is.rooted(phy)) {
    stop("Input tree must be rooted.")
  }
  if (!is.binary(phy)) {
    stop("Input tree must be binary.")
  }
  if (delta > sum(phy$edge.length)) {
    # FIXME: diameter (longest path) would be better
    stop("Input delta is greater than total tree length.")
  }
  
  root <- Ntip(phy) + 1
  step <- 0
  while (TRUE) {
    if (step == 0) {
      step <- runif(1, 0, delta)
    }
    
    child <- sample(phy$edge[which(phy$edge[,1]==root),2], 1)
    clen <- phy$edge.length[child]
    if (step > clen) {
      # step exceeds branch length
      root <- child
      step <- step - clen
      
      if (!is.element(child, phy$edge[,1])) {
        # child is terminal branch - reflect step for symmetric proposal
        phy <- reroot(phy, child)   # reroot at tip
      } 
    } 
    else {
      return(reroot(phy, child, position=clen-step))
    }
  }
}
