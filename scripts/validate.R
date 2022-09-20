# This script contains functions that we used to validate bayroot 
# on trees simulated from a compartmental model of cellular dynamics
# in the R package twt (see scripts/simulate-testdata.R)

require(chemCal)  # for inverse.predict
require(ape)
setwd("~/git/bayroot")
source("bayroot.R")
require(lubridate)

#' Fit root-to-tip regression to tree
#' @param phy:  an object of class 'phylo' (ape)
root2tip <- function(phy) {
  phy <- multi2di(phy)
  
  tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
  censored <- ifelse(grepl("Active", names(tip.dates)), tip.dates, NA)
  names(censored) <- names(tip.dates)
  
  rooted <- rtt(t=phy, tip.dates=censored)
  
  div <- node.depth.edgelength(rooted)[1:Ntip(rooted)]
  names(div) <- rooted$tip.label
  
  fit <- lm(div ~ censored)
  
  pred <- sapply(div[is.na(censored)], function(y) {
    inverse.predict(fit, y)$Prediction
  })
  
  result <- list(phy=phy, tip.dates=tip.dates, censored=censored,
                 rooted=rooted, div=div, fit=fit, pred=pred)
  class(result) <- "root2tip"
  result
}


#' Retrieve integration dates by matching tip labels in the tree to the 
#' CSV file
#' @param obj:  object of class 'root2tip'
#' @param csvfile:  path to CSV file containing actual integration dates
get.true.values <- function(obj, csvfile) {
  if (class(obj) != "root2tip") {
    stop("ERROR: expecting obj of class 'root2tip'")
  }
  # load true values (recorded as time before most recent sample)
  int.times <- read.csv(csvfile, row.names=1)
  idx <- match(row.names(int.times), 
               gsub("^(.+)_[0-9]+$", "\\1", names(obj$tip.dates)))
  true.dates <- obj$tip.dates[idx] - int.times$int.times
  
  # re-order the true values to match tip labels in the root2tip object
  idx <- match(names(obj$pred), names(true.dates))
  true.dates[idx]
}


plot.root2tip <- function(obj, true.vals, col='red', ...) {
  plot.default(obj$tip.dates, obj$div, 
       col=ifelse(is.na(obj$censored), col, 'black'), ...)
  
  abline(obj$fit)
  col.alpha <- add.alpha(col, 0.5)
  x <- obj$pred
  y <- obj$div[is.na(obj$censored)]
  
  points(x, y, col=col, pch=19, cex=0.8)
  
  segments(x0=obj$pred, x1=obj$tip.dates[is.na(obj$censored)], 
           y0=y, col=col.alpha)
  # show true dates
  points(true.vals, y, pch=3, cex=0.8, lwd=2)  
}


#' @param treefile: path to file containing Newick tree string
#' @param csvfile: path to file with 
#' @param settings:
#' @param max.date: upper limit on predicted integration dates, usually
#'                  the start of ART
fit.bayroot <- function(treefile, csvfile, settings,
                        max.date=as.Date("2000-11-01"),  # ART at 10 mo
                        #max.date=as.Date("2001-04-01"),
                        nstep=1e4, skip=10, echo=F) {
  phy <- read.tree(treefile)
  
  # modify tip labels so they can be parsed as dates
  tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
  tip.dates <- as.Date("2000-01-01") + months(tip.dates)
  
  phy$tip.label <- paste(phy$tip.label, tip.dates, sep="_")
  tip.dates <- as.Date(sapply(phy$tip.label, function(x) strsplit(x, "_")[[1]][4]))
  censored <- as.Date(ifelse(grepl("Active", phy$tip.label), tip.dates, NA), 
                      origin='1970-01-01')
  
  phy <- reroot(phy, which.min(censored))
  settings$censored <- phy$tip.label[grepl("^Latent", phy$tip.label)]
  params <- list(phy=phy, rate=0.1, origin=min(censored, na.rm=T)-1)
  
  # 1000 samples
  chain <- bayroot(nstep=nstep, skip=skip, params=params, settings=settings, echo=echo)
  
  # 200 samples
  pred.dates <- predict(chain, settings, max.date=max.date, 
                        burnin=100, thin=200)
  
  return(list(phy=phy, pred.dates=pred.dates, tip.dates=tip.dates, 
              censored=censored, chain=chain))
}


#' Restore bayroot object from logfiles
load.bayroot <- function(treefile, csvfile, settings, logfile, treelogfile,
                         max.date=as.Date('2000-11-01')) {
  phy <- read.tree(treefile)
  
  # modify tip labels so they can be parsed as dates
  tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
  tip.dates <- as.Date("2000-01-01") + months(tip.dates)
  
  phy$tip.label <- paste(phy$tip.label, tip.dates, sep="_")
  tip.dates <- as.Date(sapply(phy$tip.label, function(x) strsplit(x, "_")[[1]][4]))
  censored <- as.Date(ifelse(grepl("Active", phy$tip.label), tip.dates, NA), 
                      origin='1970-01-01')
  
  phy <- reroot(phy, which.min(censored))
  settings$censored <- phy$tip.label[grepl("^Latent", phy$tip.label)]
  
  chain <- list()
  chain$log <- read.csv(logfile)
  chain$treelog <- as.character(read.table(treelogfile, header=FALSE)$V1)
  class(chain) <- 'bayroot'
  pred.dates <- predict(chain, settings, max.date=max.date, burnin=100, thin=200)
  
  return(list(phy=phy, pred.dates=pred.dates, tip.dates=tip.dates, 
              censored=censored, chain=chain))
}


dt2months <- function(dt, refdate="2000-01-01") {
  interval(as.Date(refdate), as.Date(dt, origin="1970-01-01")) / months(1)
}


# calculate means
get.estimates <- function(obj, rt) {
  x <- sapply(obj$pred.dates, function(x) 
    quantile(x$int.date, c(0.025, 0.5, 0.975)))
  est <- apply(x, 1, dt2months)
  row.names(est) <- names(obj$pred.dates)
  
  # match to rt results
  idx <- match(
    gsub("^(.+)_[0-9]+-[0-9]+-[0-9]+$", "\\1", names(obj$pred.dates)),
    names(rt$div)
    )
  
  data.frame(est=est[,2], lo95=est[,1], hi95=est[,3], div=rt$div[idx])
}

rmse <- function(obj, true.vals) {
  sapply(1:length(true.vals), function(i) {
    pred <- dt2months(obj$pred.dates[[i]]$int.date)
    sqrt(mean((pred-true.vals[i])^2))
  })
}

