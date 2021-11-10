require(ape)
require(phytools)
require(msm)
require(ggfree)

#' shift.root
#' Reroot the input tree by moving along the branches by step size
#' drawn from uniform (0, delta).  If this step crosses an internal
#' node, go up one of the child branches at random.
#' @param phy {ape::phylo}:  input rooted tree
#' @param delta {double}:  maximum step size, s.t. proposal is U(-delta, delta)
#' 
#' @example 
#' set.seed(1)
#' phy <- rtree(5)
#' par(mfrow=c(1,2))
#' plot(phy); plot(proposal(phy, delta=0.5))
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


get.dates <- function(phy, censor=TRUE) {
  tip.dates <- as.Date(sapply(phy$tip.label, function(x) strsplit(x, "_")[[1]][4]))
  if (censor) {
    tip.dates[grepl("_DNA_", phy$tip.label)] <- NA
  }
  return (tip.dates)
}


#' lf - likelihood function
#' Calculate likelihood of divergences given input tree and sampling 
#' times.  Assume a linear relationship between divergence and time with 
#' slope {rate} and x-intercept {origin}.
#' Employ Poisson distributions centered on this trend line 
#' with variance equal to mean.
#' 
#' TODO: relax linear trend assumption with a Poisson correction (as per Jukes-
#' Cantor).  
#' TODO: relax Poisson assumption, use negative binomial distribution?
#' 
#' @param phy {ape:phylo}:  input rooted tree
#' @param origin {double}:  date at root
#' @param rate {double}:  mutation rate, i.e., slope of linear regression
#' @return {double}:  log-likelihood
lf <- function(phy, origin, rate) {
  # extract tip distances from root
  div <- node.depth.edgelength(phy)[1:Ntip(phy)]  # indexed as in phy$tip.label
  
  # time differences from origin, in days
  tip.dates <- get.dates(phy)
  delta.t <- as.integer(tip.dates - origin)
  
  # compute Poisson model
  sum(div*log(rate*delta.t) - (rate*delta.t) - lgamma(div+1), na.rm=T)
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
  max(-1e50, dunif(as.integer(origin), min=as.integer(hyper$mindate),
        max=as.integer(hyper$maxdate), log=T)) + 
    dlnorm(rate, meanlog=hyper$meanlog, sdlog=hyper$sdlog, log=T)
}


#' mh - Metropolis-Hastings sampler
#' 
#' TODO: let user specify proposal distribution parameters
#' TODO: pass file paths to write logs
#' 
#' @param nstep {integer}:  number of steps in chain sample
#' @param phy {ape:phylo}:  starting tree, rooted
#' @param tip.dates {Date}:  vector of Date objects corresponding to tip.labels
#' @param init.p {list}:  initial parameter settings
#' @param settings {list}:  hyperparameters for prior distributions and proposal settings
#' @param log.skip {integer}:  number of steps between log entries
#' @param treelog.skip {integer}:  number of steps between treelog entries
#' 
#' @return {data.frame, character}:  log, treelog
mh <- function(nstep, params, settings, log.skip=10, treelog.skip=10) {
  # deep copy
  next.params <- list(origin=params$origin, rate=params$rate, phy=params$phy)
  
  # posterior probability of initial state
  llk <- lf(params$phy, origin=params$origin, rate=params$rate)
  lprior <- prior(origin=params$origin, rate=params$rate, hyper=settings)
  lpost <- llk + lprior
  
  tip.dates <- get.dates(params$phy)
  min.date <- min(tip.dates, na.rm=T)
  
  # prepare logs
  log <- data.frame(step=0, posterior=lpost, logL=llk, prior=lprior, 
                    origin=params$origin, rate=params$rate)
  treelog <- c(write.tree(params$phy))
    
  # propagate chain sample
  for (i in 1:nstep) {
    # proposal
    next.params$phy <- shift.root(params$phy, delta=settings$root.delta)
    next.params$origin <- as.Date(round(
        rtnorm(1, mean=as.integer(params$origin), sd=settings$date.sd, 
               upper=as.integer(min.date)-1)
        ), origin='1970-01-01')
    
    rate.delta <- runif(1, min=-settings$rate.delta, max=settings$rate.delta)
    if (rate.delta <= -params$rate) {
      next.params$rate <- -(rate.delta - params$rate)  # reflect
    } else {
      next.params$rate <- params$rate + rate.delta
    }
    
    next.llk <- lf(next.params$phy, origin=next.params$origin, rate=next.params$rate)
    next.lprior <- prior(origin=next.params$origin, rate=next.params$rate, hyper=settings)
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
    }
    if (i %% treelog.skip == 0) {
      treelog <- c(treelog, write.tree(params$phy))
    }
  }
  result <- list(log=log, treelog=treelog)
  class(result) <- 'bayroot'
  return(result)
}



#' generic S3 plot for bayroot class
plot.bayroot <- function(obj, step=NA, burnin=1) {
  if (is.na(step)) {
    orig.par <- par(mfrow=c(3,2))
    end <- nrow(obj$log)
    x <- obj$log$step[burnin:end]
    
    y <- obj$log$posterior[burnin:end]
    plot(x, y, type='l', xlab='Step', ylab='Posterior')
    hist(y, main='Posterior')
    
    y <- obj$log$origin[burnin:end]
    plot(x, y, type='l', xlab='Step', ylab='Origin')
    hist(y, breaks='week', main='Origin (x-intercept)')
    
    y <- obj$log$rate[burnin:end]
    plot(x, y, type='l', xlab='Step', ylab='Rate')
    hist(y, main='Rate (slope)')
    
    par(orig.par)
  }
  else{
    phy <- read.tree(text=obj$treelog[step])
    div <- node.depth.edgelength(phy)[1:Ntip(phy)]
    tip.dates <- get.dates(phy, censor=FALSE)
    origin <- obj$log$origin[step]
    rate <- obj$log$rate[step]
    
    orig.par <- par(mfrow=c(1,2))
    plot(tree.layout(phy), mar=c(1,1,1,5), cex=0.5)
    par(mar=c(5,5,1,1))
    plot(tip.dates, div, col=ifelse(grepl("_DNA_", phy$tip.label), 'red', 'black'),
         ylim=c(0, max(div)), xlab='Sample collection date', ylab='Divergence')
    segments(x0=origin, y0=0, x1=max(tip.dates), y1= rate * (max(tip.dates) - origin))
    par(orig.par)
  }
}



#' short hand for lower incomplete gamma function, i.e.,
#' \int_0^x t^(a-1) exp(-t) dt
inc.gamma <- function(a, x) {
  pgamma(x, a) * gamma(a)
}


#' generic S3 predict for class bayroot
#' 
#' Extract sample of parameters (tree, origin, rate) from chain sample.
#' Given origin and rate, sample integration dates from the posterior probability
#' determined by sequence divergence of censored tips (DNA) for each tree.
#' 
predict.bayroot <- function(obj) {
  step <- 10  # work in progress, eventually do a sample of states

  # process tree at this step
  phy <- read.tree(text=obj$treelog[step])
  div <- node.depth.edgelength(phy)[1:Ntip(phy)]

  tip.dates <- get.dates(phy, censor=TRUE)
  min.date <- min(tip.dates, na.rm=T)
  # FIXME: actually this should be limited by sample date of censored tip
  max.date <- max(tip.dates, na.rm=T)

  origin <- obj$log$origin[step]
  dt <- max.date - min.date
  rate <- obj$log$rate[step]
  
  # calculate expected time 
  exp.t <- min.date + 1/rate * inc.gamma(div+2, rate*dt) / inc.gamma(div+1, rate*dt)
  
  list(y=origin + div/rate, x=get.dates(phy, censor=FALSE))
}


# work through test case ZM1044M - https://doi.org/10.1371/journal.ppat.1008378

setwd('~/git/bayroot')
# screened for hypermutation, aligned (MAFFT) and reconstructed tree (IQTREE)
phy <- read.tree('data/ZM1044M.fa.hyp.treefile')
#phy <- midpoint.root(phy)

# parse tip dates
tip.dates <- get.dates(phy)
phy <- reroot(phy, which.min(tip.dates))


##############  test PDF  ##############
pdf <- function(t, y, rate, t0, tmax) {
  # probability of integration time (t) given divergence (y)
  # note uniform prior on t cancels out
  L <- as.double(rate*(t-t0))
  rate * L^y * exp(-L) / inc.gamma(y+1, as.double(rate*(tmax-t0)))
}

origin <- as.Date("2006-01-26")  # sero midpoint
tip.dates <- get.dates(phy, censor=TRUE)
max.date <- max(tip.dates, na.rm=T)
x <- seq(origin, max.date, length.out=100)

y1 <- pdf(x, 0.0, 1e-4, min.date, max.date)
y2 <- pdf(x, 0.01, 1e-4, min.date, max.date)
y3 <- pdf(x, 1, 1e-4, min.date, max.date)

h <- 1/as.integer(max.date - min.date)

par(mfrow=c(1,3), cex.lab=1.2)
plot(x, y1, type='l', xlab='Sampling date', ylab='Probability density',
     main=0.001, ylim=c(0, 6.5e-4))
abline(h=h, lty=2)
plot(x, y2, type='l', col='red', xlab='Sampling date', 
     ylab='Probability density', main=0.01, ylim=c(0, 6.5e-4))
abline(h=h, lty=2)
plot(x, y3, type='l', col='blue', xlab='Sampling date', 
     ylab='Probability density', main=0.1, ylim=c(0, 6.5e-4))
abline(h=h, lty=2)



cdf <- function(t, y, rate, t0, tmax) {
  # cumulative distribution function, integrate from 0 to t
  L <- as.double(rate*(t-t0))
  inc.gamma(y+1, L) / inc.gamma(y+1, as.double(rate*(tmax-t0)))
}

y <- cdf(x, 0.01, 1e-5, min.date, max.date)

# use date of seroconversion to inform prior
settings <- list(
  # hyperparameters
  mindate=as.Date("2005-11-29"),  # last HIV -ve
  maxdate=as.Date("2006-03-25"),  # first HIV +ve
  meanlog=-10.14,  # = log(0.0144 sub/nt/yr / 365), Alizon and Fraser
  sdlog=1,  # Alizon and Fraser, 95% CI: 0.00166 - 0.0525
  
  # proposal function parameters
  root.delta=0.001,
  date.sd=10,  # days
  rate.delta=1e-6
)
# > x <- rlnorm(1e6, -10.14, 1)
# > quantile(365*x, c(0.025, 0.5, 0.975))
# 2.5%         50%       97.5% 
# 0.002027529 0.014459692 0.103403591 


init.p <- list(phy=phy, rate=1e-5, origin=min(tip.dates, na.rm=T)-1)


set.seed(2)
#results <- mh(1e3, params=init.p, settings=settings)
results <- mh(1e5, params=init.p, settings=settings, log.skip=100, treelog.skip=100)






