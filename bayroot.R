require(ape)
require(phytools)
require(msm)
require(ggfree)

#' .shift.root
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
.shift.root <- function(phy, delta=NA) {
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


#' Utility function to parse dates from sequence labels
#' @param phy {ape:phylo}:  input tree with labelled tips
#' @param delimiter {chr}:  character separating tokens in label
#' @param pos {int}:  1-index of token representing date; 
#'                    -1 indicates last token (default)
#' @param format {chr}:  AIX-style date format string; defaults to ISO format "%Y-%m-%d"
#' 
#' @export
get.dates <- function(phy, delimiter='_', pos=-1, format='%Y-%m-%d') {
  dt <- sapply(phy$tip.label, function(x) {
    tokens <- strsplit(x, delimiter)[[1]]
    if (pos == -1) {
      return(tokens[length(tokens)])
    }
    else {
      return(tokens[pos])
    }
  })
  as.Date(dt, format=format)
}


#' .lf - likelihood function
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
.lf <- function(phy, origin, rate, format="%Y-%m-%d") {
  # extract tip distances from root, indexed as in phy$tip.label
  div <- node.depth.edgelength(phy)[1:Ntip(phy)]
  
  # time differences from origin, in days
  tip.dates <- get.dates(phy, format=format)
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
.prior <- function(origin, rate, hyper) {
  #dnorm(as.integer(origin), as.integer(hyper['mean']), hyper['sd']) *
  max(-1e50, dunif(as.integer(origin), min=as.integer(hyper$mindate),
        max=as.integer(hyper$maxdate), log=T)) + 
    dlnorm(rate, meanlog=hyper$meanlog, sdlog=hyper$sdlog, log=T)
}


#' bayroot - Metropolis-Hastings sampler
#' 
#' origin = date of most recent common ancestor (root of tree), i.e.,
#'          x-intercept in root-to-tip regression
#' rate = molecular clock, expected number of substitutions 
#' TODO: pass file paths to write logs?
#' 
#' @param nstep {integer}:  number of steps in chain sample
#' @param phy {ape:phylo}:  starting tree, rooted
#' @param tip.dates {Date}:  vector of Date objects corresponding to tip.labels
#' @param init.p {list}:  initial parameter settings
#' @param settings {list}:  hyperparameters for prior distributions and proposal settings
#' @param skip {integer}:  number of steps between log entries
#' @param echo {logical}:  if TRUE, print log messages to console (stderr)
#' 
#' @return {list}:  {data frame} log, chain sample posterior, likelihood, prior,
#'                  and model parameters (origin, rate)
#'                  {character} treelog, Newick serializations of rooted trees
#'                  in chain sample.
#' @export
bayroot <- function(nstep, params, settings, skip=10, echo=FALSE) {
  # deep copy
  next.params <- list(origin=params$origin, rate=params$rate, phy=params$phy)
  
  # posterior probability of initial state
  llk <- .lf(params$phy, origin=params$origin, rate=params$rate)
  lprior <- .prior(origin=params$origin, rate=params$rate, hyper=settings)
  lpost <- llk + lprior
  
  tip.dates <- get.dates(params$phy)
  # origin cannot be more recent than first sample date
  min.date <- min(tip.dates, na.rm=T)
  
  # prepare logs
  log <- data.frame(step=0, posterior=lpost, logL=llk, prior=lprior, 
                    origin=params$origin, rate=params$rate)
  treelog <- c(write.tree(params$phy))
    
  # propagate chain sample
  for (i in 1:nstep) {
    # proposal
    next.params$phy <- .shift.root(params$phy, delta=settings$root.delta)
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
    
    next.llk <- .lf(next.params$phy, origin=next.params$origin, rate=next.params$rate)
    next.lprior <- .prior(origin=next.params$origin, rate=next.params$rate, hyper=settings)
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
    if (i %% skip == 0) {
      log <- rbind(log, list(i, lpost, llk, lprior, params$origin, params$rate))
      treelog <- c(treelog, write.tree(params$phy))
      if (echo) {
        message(i, lpost, llk, lprior, params$origin, params$rate)
      }
    }
  }
  result <- list(log=log, treelog=treelog)
  class(result) <- 'bayroot'
  return(result)
}



#' generic S3 plot for bayroot class
#' @param obj {S3}:  object of class ape::phylo
#' @param step {int}: if specified, display MCMC state at given step; otherwise
#'                    plot traces of posterior and model parameters.
#' @param burnin {int}: number of steps in chain sample to discard as burnin.
#' @export
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
    tip.dates <- get.dates(phy)
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
#' @param a {double}: gamma exponent parameter
#' @param x {double}: upper limit of integration
.inc.gamma <- function(a, x, log=FALSE) {
  if (log) {
    pgamma(x, a, log=TRUE) + lgamma(a)
  } else {
    pgamma(x, a) * gamma(a)  
  }
}


#' probability distribution function for integration time (t)
#' given divergence (y)
.pdfunc <- function(t, y, rate, t0, tmax) {
  L <- rate*(as.integer(t)-as.integer(t0))
  exp(log(rate) + y*log(L) - L - 
        .inc.gamma(y+1, rate*(as.integer(tmax)-as.integer(t0)), log=T))
}

#' generate random deviates by rejection sampling
.sample.pdfunc <- function(y, rate, t0, tmax, max.tries=1e3) {
  # calculate maximum value for PDF
  # TODO: is there a closed form solution?
  f <- function(t) { -.pdfunc(t, y, rate, t0, tmax) }  # neg for minimization
  mid.point <- t0 + (tmax-t0)/2
  res <- optim(mid.point, f, method='Brent', lower=t0, upper=tmax)
  fmax <- -res$value
  
  # rejection sampling
  f <- function(t) { .pdfunc(t, y, rate, t0, tmax) }
  tries <- 0
  while(TRUE) {
    tries <- tries + 1
    # sample random date
    s <- as.Date(runif(1, min=t0, max=tmax), origin='1970-01-01')
    fs <- f(s)  # evaluate pdf at t=s
    if (runif(1, 0, fmax) < fs) {
      break  # accept
    }
    if (tries > max.tries) {
      warning("Exceeded maximum attempts", tries, "in .sample.pdfunc()")
      return(NA)
    }
  }
  return(s)
}


#' generic S3 predict for class bayroot
#' 
#' Extract sample of parameters (tree, origin, rate) from chain sample.
#' Given origin and rate, sample integration dates from the posterior probability
#' determined by sequence divergence of censored tips (DNA) for each tree.
#' 
#' @param obj {ape:phylo}:  S3 object of class ape:phylo
#' @param censored {character}:  tip labels to predict dates for
#' @param max.date {Date}:  set upper bound on date estimates for censored tips.
#'                          Defaults to most recent date for uncensored tips.
#' @param burning {integer}:  number of steps to discard from start of chain sample
#' @param thin {integer}:  number of steps to sample at regular intervals from 
#'                         post-burning chain
#' @param delimiter {character}:  character separating tokens in label
#' @param pos {integer}:  1-index of token representing date; -1 indicates last token (default)
#' @param format {character}:  AIX-style date format string; defaults to ISO format "%Y-%m-%d"
#' @return {data.frame}:  sampled dates as {double} values; each row corresponds
#'                        to a step, and each column to a censored tip.
#' @export
predict.bayroot <- function(obj, censored, max.date=NA, burnin=10, thin=100,
                            delimiter="_", pos=-1, format="%Y-%m-%d") {
  rows <- as.integer(seq(burnin, nrow(obj$log), length.out=thin))
  
  phy <- read.tree(text=obj$treelog[1])
  tip.dates <- get.dates(phy, delimiter=delimiter, pos=pos, format=format)
  if (is.na(max.date)) {
    # limit to uncensored tips
    max.date <- max(tip.dates[!is.element(phy$tip.label, censored)], na.rm=T)
  }
  
  # prepare output container
  res <- matrix(NA, nrow=length(rows), ncol=length(censored))
  res <- as.data.frame(res)
  names(res) <- censored
  row.names(res) <- rows
  
  for (i in 1:length(rows)) {
    step <- rows[i]
    
    # unpack state
    origin <- obj$log$origin[step]
    rate <- obj$log$rate[step]
    phy <- read.tree(text=obj$treelog[step])
    
    # tip label order changes with re-rooting and calling read.tree()
    labels <- phy$tip.label
    dates <- get.dates(phy, delimiter=delimiter, pos=pos, format=format)
    
    # sample integration times for censored tips given divergence
    div <- node.depth.edgelength(phy)[1:Ntip(phy)]
    samp <- sapply(which(is.element(labels, censored)), function(i) {
      # in case censored tip sampled before last uncensored tip (max.date)
      this.max.date <- min(max.date, dates[i])
      .sample.pdfunc(div[i], rate, origin, this.max.date)
      })
    names(samp) <- labels[is.element(labels, censored)]
    
    # append sampled dates to container
    res[i, ] <- samp[match(names(samp), names(res))]
  }
  return(res)
}



