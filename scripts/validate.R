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
  mindate=as.Date("1999-12-01"), 
  maxdate=as.Date("2000-04-01"),  # first sample is 3 months after origin
  #maxdate=as.Date("2000-12-01"),  # 11 mo after origin
  
  meanlog=-5, sdlog=2,  # rate
  
  # proposal parameters
  root.delta=0.01, #2*mean(phy$edge.length),
  date.sd=10,  # days, origin proposal
  rate.delta=0.01
)


fit.bayroot <- function(treefile, csvfile, settings,
                        max.date=as.Date("2000-11-01"),  # ART at 10 mo
                        #max.date=as.Date("2001-04-01"),
                        nstep=1e4, skip=10) {
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
  chain <- bayroot(nstep=nstep, skip=skip, params=params, settings=settings)
  
  # 200 samples
  pred.dates <- predict(chain, settings, max.date=max.date, 
                        burnin=100, thin=200)
  
  return(list(phy=phy, pred.dates=pred.dates, tip.dates=tip.dates, 
              censored=censored, chain=chain))
}


load.bayroot <- function(treefile, csvfile, settings, logfile, treelogfile) {
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
  chain$treelog <- as.character(read.table(treelogfile, header=FALSE))
  class(chain) <- 'bayroot'
  pred.dates <- predict(chain, settings, max.date=max.date, burnin=100, thin=200)
  
  return(list(phy=phy, pred.dates=pred.dates, tip.dates=tip.dates, 
              censored=censored, chain=chain))
}


dt2months <- function(dt, refdate="2000-01-01") {
  interval(as.Date(refdate), as.Date(dt, origin="1970-01-01")) / months(1)
}
#plot(rt, true.vals, ylim=c(0, 0.12))
#for (temp in res$pred.dates) {
#  points(x=dt2months(temp$int.date), y=temp$div/settings$seq.len, pch=19, cex=0.5, 
#         col=rgb(0,0,1,0.2))  
#}

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


files <- Sys.glob("data/latent2.*.ft2.nwk")


if (FALSE) {
  # try one replicate first
  tf <- files[47]
  cf <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".times.csv", tf)
  phy <- read.tree(tf)
  rt <- root2tip(phy)
  true.vals <- get.true.values(rt, cf)
  
  res <- fit.bayroot(tf, cf, settings=settings, nstep=2e4, skip=20)
  
  # display results
  #pdf(file="~/slides/img/latent2.1.compare.pdf", width=5, height=4)
  par(mar=c(5,5,1,1))
  plot(rt, true.vals=true.vals, xlab="Collection date (months since origin)",
       ylab="Divergence", xlim=c(0, 20), 
       #ylim=c(0.0025, 0.075), 
       las=1, cex.axis=0.8)
  #text(rt$tip.dates[is.na(rt$censored)], rt$div[is.na(rt$censored)], 
  #     label=names(rt$tip.dates)[is.na(rt$censored)], cex=0.5, adj=1)
  
  est <- get.estimates(res, rt)
  points(est$est, est$div, pch=19, col=rgb(0,0,1,0.5), cex=0.8)
  
  segments(x0=est$lo95, x1=est$hi95, y0=est$div, col=rgb(0,0,1,0.3), lwd=5)
  abline(v=10, lty=2)
  #abline(v=15, lty=2)
  legend(x=1, y=0.035, legend=c("RTT", "bayroot", "true date"), cex=0.8,
         col=c('red', 'blue', 'black'), pch=c(19, 19, 3), pt.lwd=2, bty='n')
  #dev.off()
}


# batch processing
results <- data.frame(filename=files, rtt=NA, bayroot=NA)
for(i in 1:length(files)) {
  tf <- files[i]
  print(tf)
  cf <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".times.csv", tf)
  logfile <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".log.csv", tf)
  treefile <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".treelog.nwk", tf)
  
  phy <- read.tree(tf)
  
  # root to tip
  rt <- root2tip(phy)
  true.vals <- get.true.values(rt, cf)
  results$rtt[i] <- sqrt(mean((rt$pred-true.vals)^2))
  
  # Bayesian
  res <- fit.bayroot(tf, cf, settings=settings, nstep=2e4, skip=20)
  write.csv(res$chain$log, file=logfile)
  con <- file(treefile)
  writeLines(res$chain$treelog, con)
  close(con)
  #res <- read.csv(logfile)
  
  est <- get.estimates(res, rt)
  results$bayroot[i] <- sqrt(mean((est$est - true.vals)^2))  
}

write.csv(results, "~/git/bayroot/sim-rtt-results2.csv")

# concordance in RMSE
par(mar=c(5,5,1,1))
plot(results$rtt, results$bayroot, xlim=c(0, 4), ylim=c(0, 4))
abline(a=0, b=1, lty=2)

wilcox.test(results$rtt, results$bayroot, paired=T)
#Wilcoxon signed rank test with continuity correction
#
#data:  results$rtt and results$bayroot
#V = 1008, p-value = 0.0003547
#alternative hypothesis: true location shift is not equal to 0

png("~/slides/img/bayroot-slopegraph2.png", width=4*150, height=6*150, res=150)
par(family='Palatino', mar=c(4,4,1,1))
slopegraph(results$rtt, results$bayroot, names.arg=c('RTT', 'bayroot'), 
           colorize=T, type='l', ylab="Root mean square error", shim=0)
dev.off()

