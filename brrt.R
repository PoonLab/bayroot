require(ape)
require(phytools)
require(msm)


#' shift.root
#' Reroot the input tree by moving along the branches by step size
#' drawn from uniform (0, delta).  If this step crosses an internal
#' node, go up the left branch with probability 0.5.
#' @param phy {ape::phylo}:  input rooted tree
#' @param delta {double}:  maximum step size, s.t. proposal is U(-delta, delta)
shift.root <- function(phy, delta=NA) {
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
    #print(paste(root, child, clen, step))
    
    if (step > clen) {
      # step exceeds branch length
      root <- child
      step <- step - clen
    
      if (!is.element(child, phy$edge[,1])) {
        # child is terminal branch - reflect step for symmetric proposal
        phy <- reroot(phy, child)   # reroot at tip
        
        # update child
        prev.root <- phy$edge[,1][phy$edge[,2] == child]
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
#phy <- read.tree(text="(t2:0.06178627047,((t1:0.6870228467,(t3:0.76984142,t4:0.4976992421):0.3841037182):0.1765567525,t5:0.7176185083):0.2059745749);")
#par(mfrow=c(1,2))
#plot(phy); plot(proposal(phy, delta=0.5))


# TODO: Metropolis-Hastings - from another package?

#' lf - likelihood function
#' Calculate likelihood of divergences given input tree and sampling 
#' times.  Assume a linear relationship between divergence and time.
#' TODO: relax this assumption with a Poisson correction (as per Jukes-
#' Cantor).  Employ Poisson distributions centered on this trend line 
#' with variance equal to mean.
#' TODO: relax this assumption, use negative binomial distribution?
#' @param phy {ape:phylo}:  input rooted tree
#' @param origin {double}:  date at root
#' @param rate {double}:  mutation rate, i.e., slope of linear regression
#' @param tip.dates {numeric}:  
lf <- function(phy, origin, rate, tip.dates) {
  # extract tip distances from root
  div <- node.depth.edgelength(phy)[1:Ntip(phy)]  # indexed as in phy$tip.label
  
  # time differences from origin, in days
  delta.t <- as.integer(tip.dates - origin)
  
  # compute Poisson model
  sum(div * log(rate*delta.t) - (rate*delta.t) - lgamma(div*1), na.rm=T)
}

#' prior probability
#' We assume origin and rate are independent.
#' @param origin {double}:  origin date
#' @param rate {double}:  mutation rate
#' @param hyper {list}:  hyperparameters for prior distributions
#'    origin:  dnorm(mean, sd)
#'    rate:  dlnorm(meanlog, sdlog)
prior <- function(origin, rate, hyper) {
  dnorm(as.integer(origin), as.integer(hyper['mean']), hyper['sd']) *
    dlnorm(hyper['meanlog'], hyper['sdlog'])
}

# use date of seroconversion to inform prior
hyper <- list(
  mean=min(tip.dates, na.rm=T), 
  sd=30,  # days
  meanlog=-10,
  sdlog=2
  )


init.p <- list(rate=1e-5, origin=min(tip.dates, na.rm=T)-1)



#' mh - Metropolis-Hastings sampler
#' TODO: let user specify proposal distribution parameters
#' @param nstep {integer}:  number of steps in chain sample
#' @param phy {ape:phylo}:  starting tree, rooted
#' @param tip.dates {Date}:  vector of Date objects corresponding to tip.labels
#' @param init.p {list}:  initial parameter settings
#' @param hyper {list}:  hyperparameters for prior distributions
mh <- function(nstep, phy, tip.dates, init.p, hyper) {
  # unpack parameters
  params <- list(origin=init.p$origin, rate=init.p$rate, phy=phy)
  next.params <- list(origin=init.p$origin, rate=init.p$rate, phy=phy)
  
  # posterior probability of initial state
  pp <- lf(params$phy, origin=params$origin, rate=params$rate, 
           tip.dates=tip.dates)
  min.date <- min(tip.dates, na.rm=T)
  
  # propagate chain sample
  for (i in 1:nstep) {
    print(paste(i, pp, params$origin, params$rate))
    
    # proposal
    u <- runif(1)
    if (u < 0.5) {
      next.params$phy <- shift.root(params$phy, delta=0.001)
    } else if (u < 0.75) {
      next.params$origin <- as.Date(round(
        rtnorm(1, mean=as.integer(params$origin), sd=7, upper=as.integer(min.date)-1)
        ), origin='1970-01-01')
    } else {
      next.params$rate <- rlnorm(1, meanlog=log(params$rate), sdlog=0.1)
    }
    
    next.pp <- lf(next.params$phy, origin=next.params$origin, 
                  rate=next.params$rate, tip.dates=tip.dates)
    #print(paste("next.pp", next.pp))
    if (is.infinite(next.pp)) {
      stop(next.params)
    }
    
    ratio <- next.pp/pp
    if (ratio >= 1 | runif(1) < ratio) {
      # accept proposal
      params <- next.params
      pp <- next.pp
    }
  }
}



setwd('~/git/brrt')
phy <- read.tree('data/ZM1044M.fa.hyp.treefile')
phy <- midpoint.root(phy)

# parse tip dates
tip.dates <- as.Date(sapply(phy$tip.label, function(x) strsplit(x, "_")[[1]][4]))
tip.dates[grepl("_DNA_", phy$tip.label)] <- NA

set.seed(2)
mh(100, phy, tip.dates, init.p, hyper)

