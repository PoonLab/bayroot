# apply ML root-to-tip regression to simulated data
# 
require(chemCal)  # for inverse.predict

fit.rtt <- function(treefile, csvfile, plot=FALSE) {
  phy <- read.tree(treefile)
  phy <- multi2di(phy)

  # fit RTT model
  tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
  censored <- ifelse(grepl("Active", names(tip.dates)), tip.dates, NA)
  rooted <- rtt(t=phy, tip.dates=censored)
  div <- node.depth.edgelength(rooted)[1:Ntip(rooted)]
  mod <- lm(div ~ censored)
  
  # generate predictions
  pred.dates <- sapply(div[is.na(censored)], function(y) {
    inverse.predict(mod, y)$Prediction
  })

  # load true values (recorded as time before most recent sample)
  int.times <- read.csv(csvfile, row.names=1)
  idx <- match(row.names(int.times), 
               gsub("^(.+)_[0-9]+$", "\\1", names(tip.dates)))
  true.dates <- tip.dates[idx] - int.times$int.times
  idx <- match(names(true.dates), rooted$tip.label[is.na(censored)])
  true.dates <- true.dates[idx]
  
  if (plot) {
    par(mar=c(5,5,1,1))
    plot(tip.dates, div, xlim=c(0, 20), ylim=range(div),
         pch=ifelse(is.na(censored), 19, 1))
    
    # show predicted times
    abline(mod)
    points(pred.dates, div[is.na(censored)], col='red', pch=19, cex=0.8)
    segments(x0=pred.dates, x1=tip.dates[is.na(censored)], 
             y0=div[is.na(censored)], col='red')
    # show true dates
    points(true.dates, div[is.na(censored)], col='blue', pch=3, cex=0.8)
  }
  
  # root-mean square error
  sqrt(mean((pred.dates-true.dates[idx])^2))
}

setwd("~/git/bayroot")
files <- Sys.glob("data/latent1.*.ft2.nwk")

for(tf in files) {
  cf <- gsub("\\.cens\\.nwk\\.fas\\.ft2\\.nwk", ".times.csv", tf)
  print(fit.rtt(tf, cf))
}
