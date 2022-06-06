# apply ML root-to-tip regression to simulated data
# 
require(chemCal)  # for inverse.predict
require(ape)
setwd("~/git/bayroot")
source("bayroot.R")
require(lubridate)

#' Fit root-to-tip regression to tree
root2tip <- function(phy) {
  phy <- multi2di(phy)
  
  tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
  censored <- ifelse(grepl("Active", names(tip.dates)), tip.dates, NA)
  names(censored) <- names(tip.dates)
  
  rooted <- rtt(t=phy, tip.dates=censored)
  
  div <- node.depth.edgelength(rooted)[1:Ntip(rooted)]
  names(div) <- rooted$tip.label
  
  fit <- lm(div ~ censored)
  
  pred <- sapply(obj$div[is.na(obj$censored)], function(y) {
    inverse.predict(obj$fit, y)$Prediction
  })
  
  result <- list(phy=phy, tip.dates=tip.dates, censored=censored,
                 rooted=rooted, div=div, fit=fit, pred=pred)
  class(result) <- "root2tip"
  result
}


get.true.values <- function(obj, csvfile) {
  if (class(obj) != "root2tip") {
    stop("ERROR: expecting obj of class 'root2tip'")
  }
  # load true values (recorded as time before most recent sample)
  int.times <- read.csv(csvfile, row.names=1)
  idx <- match(row.names(int.times), 
               gsub("^(.+)_[0-9]+$", "\\1", names(obj$tip.dates)))
  
  true.dates <- obj$tip.dates[idx] - int.times$int.times
  idx <- match(names(obj$pred), names(true.dates))
  true.dates[idx]
}


plot.root2tip <- function(obj, true.vals, ...) {
  plot(obj$tip.dates, obj$div, xlim=c(0, 20), ylim=range(obj$div),
       col=ifelse(is.na(obj$censored), 'red', 'black'), ...)
  
  abline(obj$fit)
  red <- rgb(1,0,0,0.5)
  points(obj$pred, obj$div[is.na(obj$censored)], col=red, pch=19, cex=0.8)
  segments(x0=obj$pred, x1=obj$tip.dates[is.na(obj$censored)], 
           y0=obj$div[is.na(obj$censored)], col=red)
  # show true dates
  points(true.vals, obj$div[is.na(obj$censored)], pch=3, cex=0.8, lwd=2)  
}



fit.bayroot <- function(treefile, csvfile, nstep=1e4, skip=10) {
  phy <- read.tree(treefile)
  
  # modify tip labels so they can be parsed as dates
  tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
  tip.dates <- as.Date("2000-01-01") + months(tip.dates)
  phy$tip.label <- paste(phy$tip.label, tip.dates, sep="_")
  
  tip.dates <- as.Date(sapply(phy$tip.label, function(x) strsplit(x, "_")[[1]][4]))
  censored <- as.Date(ifelse(grepl("Active", phy$tip.label), tip.dates, NA), 
                      origin='1970-01-01')
  
  phy <- reroot(phy, which.min(censored))
  params <- list(phy=phy, rate=0.1, origin=min(censored, na.rm=T)-1)
  
  settings <- list(
    # hyperpaameters
    mindate=as.Date("2000-01-01"), maxdate=max(censored, na.rm=T),  # origin
    meanlog=0, sdlog=1,  # rate
    
    # proposal parameters
    root.delta=0.1*median(phy$edge.length),
    date.sd=10,  # days
    rate.delta=0.001
    )
  
  # 1000 samples
  chain <- bayroot(nstep=nstep, skip=skip, params=params, settings=settings)
  
  # 200 samples
  pred.dates <- predict(chain, phy$tip.label[is.na(censored)], burnin=100, thin=180)
  
  return(list(phy=phy, pred.dates=pred.dates, tip.dates=tip.dates, 
              censored=censored, log=chain$log, treelog=chain$treelog))
}




files <- Sys.glob("data/latent1.*.ft2.nwk")
tf <- files[1]
cf <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".times.csv", tf)

# try one replicate first
phy <- read.tree(tf)
rt <- root2tip(phy)
true.vals <- get.true.values(rt, cf)

rmse.rtt <- sqrt(mean((pred.vals-true.vals)^2))

res <- fit.bayroot(tf, cf)
bay.mean <- as.Date(apply(res$pred.dates, 2, mean), origin="1970-01-01")
bay.lo <- as.Date(apply(res$pred.dates, 2, function(x) quantile(x, 0.025)), 
                  origin="1970-01-01")
bay.hi <- as.Date(apply(res$pred.dates, 2, function(x) quantile(x, 0.975)), 
                  origin="1970-01-01")

# convert to number of months since 2000-01-01
bay.mean <- interval(as.Date("2000-01-01"), bay.mean) / months(1)
bay.lo <- interval(as.Date("2000-01-01"), bay.lo) / months(1)
bay.hi <- interval(as.Date("2000-01-01"), bay.hi) / months(1)

rmse.bay <- sqrt(mean((bay.mean - true.vals)^2))

names(bay.mons) <- gsub("^(.+)_[0-9]+-[0-9]+-[0-9]+$", "\\1", names(bay.dates))

plot(rt, true.vals=true.vals, xlab="Collection date (months since origin)",
     ylab="Divergence")
points(bay.mons, rt$div[is.na(rt$censored)], pch=19, col=rgb(0,0,1,0.5), cex=0.8)
segments(x0=bay.lo, x1=bay.hi, y0=rt$div[is.na(rt$censored)], col=rgb(0,0,1,0.5))


results <- data.frame(filename=files, rtt=NA, bayroot=NA)
for(i in 1:length(files)) {
  tf <- files[i]
  print(tf)
  cf <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".times.csv", tf)
  phy <- read.tree(tf)
  
  # root to tip
  rt <- root2tip(phy)
  true.vals <- get.true.values(rt, cf)
  results$rtt[i] <- sqrt(mean((pred.vals-true.vals)^2))
  
  # Bayesian
  res <- fit.bayroot(tf, cf)
  bay.mean <- as.Date(apply(res$pred.dates, 2, mean), origin="1970-01-01")
  bay.mean <- interval(as.Date("2000-01-01"), bay.mean) / months(1)
  results$bayroot[i] <- sqrt(mean((bay.mean - true.vals)^2))  
}
