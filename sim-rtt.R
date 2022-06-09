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
  
  pred <- sapply(div[is.na(censored)], function(y) {
    inverse.predict(fit, y)$Prediction
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
  plot.default(obj$tip.dates, obj$div, 
       col=ifelse(is.na(obj$censored), 'red', 'black'), ...)
  
  abline(obj$fit)
  red <- rgb(1,0,0,0.5)
  x <- obj$pred
  y <- obj$div[is.na(obj$censored)]
  
  points(x, y, col=red, pch=19, cex=0.8)
  
  segments(x0=obj$pred, x1=obj$tip.dates[is.na(obj$censored)], 
           y0=y, col=red)
  # show true dates
  points(true.vals, y, pch=3, cex=0.8, lwd=2)  
}


settings <- list(
  seq.len=1233,  # AY772699
  format="%Y-%m-%d",
  
  # hyperparameters
  mindate=as.Date("1999-01-01"), 
  maxdate=as.Date("2000-04-01"),  # origin
  
  meanlog=-5, sdlog=2,  # rate
  
  # proposal parameters
  root.delta=2*mean(phy$edge.length),
  date.sd=10,  # days, origin proposal
  rate.delta=0.01
)


fit.bayroot <- function(treefile, csvfile, settings, nstep=1e4, skip=10) {
  phy <- read.tree(treefile)
  
  # modify tip labels so they can be parsed as dates
  tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
  tip.dates <- as.Date("2000-01-01") + months(tip.dates)
  
  phy$tip.label <- paste(phy$tip.label, tip.dates, sep="_")
  settings$censored <- phy$tip.label[grepl("^Latent", phy$tip.label)]
  
  tip.dates <- as.Date(sapply(phy$tip.label, function(x) strsplit(x, "_")[[1]][4]))
  censored <- as.Date(ifelse(grepl("Active", phy$tip.label), tip.dates, NA), 
                      origin='1970-01-01')
  
  phy <- reroot(phy, which.min(censored))
  params <- list(phy=phy, rate=0.1, origin=min(censored, na.rm=T)-1)
  
  # 1000 samples
  chain <- bayroot(nstep=nstep, skip=skip, params=params, settings=settings)
  
  # 200 samples
  pred.dates <- predict(chain, settings, max.date=as.Date("2000-11-01"), 
                        burnin=100, thin=200)
  
  return(list(phy=phy, pred.dates=pred.dates, tip.dates=tip.dates, 
              censored=censored, chain=chain))
}




files <- Sys.glob("data/latent1.*.ft2.nwk")
tf <- files[2]
cf <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".times.csv", tf)

# try one replicate first
phy <- read.tree(tf)
rt <- root2tip(phy)
true.vals <- get.true.values(rt, cf)

rmse.rtt <- sqrt(mean((rt$pred-true.vals)^2))

res <- fit.bayroot(tf, cf, settings=settings, nstep=2e4, skip=20)

#temp <- res$pred.dates[[1]]
dt2months <- function(dt, refdate="2000-01-01") {
  interval(as.Date(refdate), as.Date(dt, origin="1970-01-01")) / months(1)
}
for (temp in res$pred.dates) {
  points(x=dt2months(temp$int.date), y=temp$div/settings$seq.len, pch=19, cex=0.5, 
         col=rgb(0,0,1,0.2))  
}


# calculate means
bay.mean <- as.Date(sapply(res$pred.dates, function(x) mean(x$int.date)),
                    origin="1970-01-01")
bay.lo <- as.Date(sapply(res$pred.dates, function(x) quantile(x$int.date, 0.025)),
                  origin="1970-01-01")
bay.hi <- as.Date(sapply(res$pred.dates, function(x) quantile(x$int.date, 0.975)),
                  origin="1970-01-01")

# convert to number of months since 2000-01-01
bay.mean <- interval(as.Date("2000-01-01"), bay.mean) / months(1)
bay.lo <- interval(as.Date("2000-01-01"), bay.lo) / months(1)
bay.hi <- interval(as.Date("2000-01-01"), bay.hi) / months(1)

rmse.bay <- sqrt(mean((bay.mean - true.vals)^2))

bay.div <- sapply(res$pred.dates, function(x) mean(x$div) / settings$seq.len)

# display results
#pdf(file="latent1.1.compare.pdf", width=5, height=5)
plot(rt, true.vals=true.vals, xlab="Collection date (months since origin)",
     ylab="Divergence")
points(bay.mean, bay.div, pch=19, col=rgb(0,0,1,0.5), cex=0.8)
segments(x0=bay.lo, x1=bay.hi, y0=bay.div, col=rgb(0,0,1,0.5))
abline(v=10, lty=2)
legend(x=13, y=0.1, legend=c("RTT", "bayroot", "true date"), 
       col=c('red', 'blue', 'black'), pch=c(19, 19, 3), pt.lwd=2)
#dev.off()



# batch processing
results <- data.frame(filename=files, rtt=NA, bayroot=NA)
for(i in 1:length(files)) {
  tf <- files[i]
  print(tf)
  cf <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".times.csv", tf)
  phy <- read.tree(tf)
  
  # root to tip
  rt <- root2tip(phy)
  true.vals <- get.true.values(rt, cf)
  results$rtt[i] <- sqrt(mean((rt$pred-true.vals)^2))
  
  # Bayesian
  res <- fit.bayroot(tf, cf, settings=settings, nstep=2e4, skip=20)
  bay.mean <- as.Date(sapply(res$pred.dates, function(x) mean(x$int.date)),
                      origin="1970-01-01")
  bay.mean <- interval(as.Date("2000-01-01"), bay.mean) / months(1)
  results$bayroot[i] <- sqrt(mean((bay.mean - true.vals)^2))  
}
