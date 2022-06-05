# apply ML root-to-tip regression to simulated data
# 
require(chemCal)  # for inverse.predict
require(ape)
setwd("~/git/bayroot")
source("bayroot.R")
require(lubridate)

fit.rtt <- function(treefile) {
  phy <- read.tree(treefile)
  phy <- multi2di(phy)

  # fit RTT model
  tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
  censored <- ifelse(grepl("Active", names(tip.dates)), tip.dates, NA)
  rooted <- rtt(t=phy, tip.dates=censored)
  div <- node.depth.edgelength(rooted)[1:Ntip(rooted)]
  mod <- lm(div ~ censored)
  
  result <- list(tip.dates=tip.dates, censored=censored, rooted=rooted, 
                 div=div, mod=mod, pred.dates=pred.dates)
  class(result) <- "fitrtt"
  result
}


predict.fitrtt <- function(obj, csvfile) {
  # generate predictions
  sapply(obj$div[is.na(obj$censored)], function(y) {
    inverse.predict(obj$mod, y)$Prediction
  })
}


plot.fitrtt <- function(obj) {
  plot(obj$tip.dates, obj$div, xlim=c(0, 20), ylim=range(obj$div),
       pch=ifelse(is.na(obj$censored), 19, 1))
  abline(obj$mod)

  # show predicted times
  pred.dates <- predict(obj)
  points(pred.dates, obj$div[is.na(obj$censored)], col='red', pch=19, cex=0.8)
  segments(x0=pred.dates, x1=obj$tip.dates[is.na(obj$censored)], 
           y0=obj$div[is.na(obj$censored)], col='red')
  
  # load true values (recorded as time before most recent sample)
  int.times <- read.csv(csvfile, row.names=1)
  idx <- match(row.names(int.times), 
               gsub("^(.+)_[0-9]+$", "\\1", names(obj$tip.dates)))
  true.dates <- obj$tip.dates[idx] - int.times$int.times
  idx <- match(names(true.dates), obj$rooted$tip.label[is.na(obj$censored)])
  true.dates <- true.dates[idx]
  points(true.dates, obj$div[is.na(obj$censored)], col='blue', pch=3, cex=0.8)  
  
  # root-mean square error
  sqrt(mean((pred.dates-true.dates[idx])^2))
}



fit.bayroot <- function(treefile, csvfile, nstep=1e4, skip=100) {
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
  
  chain <- mh(nstep=nstep, skip=skip, params=params, settings=settings)
  pred.dates <- predict.bayroot(chain, phy$tip.label[is.na(censored)])
  return(list(pred.dates=pred.dates, log=chain$log, treelog=chain$treelog))
}



files <- Sys.glob("data/latent1.*.ft2.nwk")

for(tf in files) {
  cf <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".times.csv", tf)
  print(fit.rtt(tf, cf))
}
