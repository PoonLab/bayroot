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
  #dnorm(as.integer(origin), as.integer(hyper['mean']), hyper['sd']) *
  max(0, dunif(as.integer(origin), min=as.integer(hyper$mindate),
        max=as.integer(hyper$maxdate), log=T)) + 
    dlnorm(rate, meanlog=hyper$meanlog, sdlog=hyper$sdlog, log=T)
}


init.p <- list(rate=1e-5, origin=min(tip.dates, na.rm=T)-1)



#' mh - Metropolis-Hastings sampler
#' 
#' TODO: let user specify proposal distribution parameters
#' TODO: pass file paths to write logs
#' 
#' @param nstep {integer}:  number of steps in chain sample
#' @param phy {ape:phylo}:  starting tree, rooted
#' @param tip.dates {Date}:  vector of Date objects corresponding to tip.labels
#' @param init.p {list}:  initial parameter settings
#' @param hyper {list}:  hyperparameters for prior distributions
#' @param log.skip {integer}:  number of steps between log entries
#' @param treelog.skip {integer}:  number of steps between treelog entries
#' 
#' @return {data.frame, character}:  log, treelog
mh <- function(nstep, phy, tip.dates, init.p, hyper, log.skip=10, treelog.skip=10) {
  # unpack parameters
  params <- list(origin=init.p$origin, rate=init.p$rate, phy=phy)
  next.params <- list(origin=init.p$origin, rate=init.p$rate, phy=phy)
  
  # posterior probability of initial state
  llk <- lf(params$phy, origin=params$origin, rate=params$rate, 
           tip.dates=tip.dates)
  lprior <- prior(origin=params$origin, rate=params$rate, hyper=hyper)
  lpost <- llk + lprior
  min.date <- min(tip.dates, na.rm=T)
  
  # prepare logs
  log <- data.frame(step=0, posterior=lpost, logL=llk, prior=lprior, 
                    origin=params$origin, rate=params$rate)
  treelog <- c(write.tree(params$phy))
    
  # propagate chain sample
  for (i in 1:nstep) {
    # proposal
    next.params$phy <- shift.root(params$phy, delta=0.001)
    next.params$origin <- as.Date(round(
        rtnorm(1, mean=as.integer(params$origin), sd=10, upper=as.integer(min.date)-1)
        ), origin='1970-01-01')
    
    rate.delta <- runif(1, min=-1e-6, max=1e-6)
    if (rate.delta <= -params$rate) {
      next.params$rate <- -(rate.delta - params$rate)  # reflect
    } else {
      next.params$rate <- params$rate + rate.delta
    }
    
    next.llk <- lf(next.params$phy, origin=next.params$origin, 
                  rate=next.params$rate, tip.dates=tip.dates)
    next.lprior <- prior(origin=next.params$origin, rate=next.params$rate, hyper=hyper)
    next.lpost <- next.llk + next.lprior
    if (is.infinite(next.lpost)) {
      stop(next.params)
    }
    
    ratio <- exp(next.lpost - lpost)
    if (ratio >= 1 | runif(1) < ratio) {
      # accept proposal
      params <- next.params
      lpost <- next.lpost
      llk <- next.llk
      lprior <- next.lprior
    }
    
    # update logs
    if (i %% log.skip == 0) {
      log <- rbind(log, list(i, lpost, llk, lprior, params$origin, params$rate))
      treelog <- c(treelog, write.tree(params$phy))
    }
  }
  return(list(log=log, treelog=treelog))
}


# work through test case ZM1044M - https://doi.org/10.1371/journal.ppat.1008378

setwd('~/git/brrt')
# screened for hypermutation, aligned (MAFFT) and reconstructed tree (IQTREE)
phy <- read.tree('data/ZM1044M.fa.hyp.treefile')
phy <- midpoint.root(phy)

# parse tip dates
tip.dates <- as.Date(sapply(phy$tip.label, function(x) strsplit(x, "_")[[1]][4]))
tip.dates[grepl("_DNA_", phy$tip.label)] <- NA

# use date of seroconversion to inform prior
hyper <- list(
  mindate=as.Date("2005-11-29"),  # last HIV -ve
  maxdate=as.Date("2006-03-25"),  # first HIV +ve
  meanlog=-10.14,  # Alizon and Fraser, 10^-1.84 sub/nt/year 
  sdlog=1 # (Alizon and Fraser, 95% CI: 10^-2.78 - 10^-1.28)
)

set.seed(2)
results <- mh(1e3, phy, tip.dates, init.p, hyper)
#results <- mh(1e5, phy, tip.dates, init.p, hyper, log.skip=100, treelog.skip=1000)

plot(results$log$step, results$log$posterior, type='l')
plot(results$log$step, results$log$origin, type='l')
plot(results$log$step, results$log$rate, type='l')
