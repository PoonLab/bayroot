source("scripts/validate.R")

# MAIN ROUTINE
files <- Sys.glob("data/latent1.*.ft2.nwk")

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
  
  # Bayesian sampling
  res <- fit.bayroot(tf, cf, settings=settings, nstep=2e4, skip=20)
  
  # write outputs to files
  write.csv(res$chain$log, file=logfile)
  con <- file(treefile)
  writeLines(res$chain$treelog, con)
  close(con)
  
  #res <- read.csv(logfile)
  est <- get.estimates(res, rt)
  results$bayroot[i] <- sqrt(mean((est$est - true.vals)^2))  
}

write.csv(results, "~/git/bayroot/sim-rtt-results2.csv")
#results <- read.csv("~/git/bayroot/working/sim-rtt-results2.csv")

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

#png("~/slides/img/bayroot-slopegraph2.png", width=4*150, height=6*150, res=150)
pdf("~/papers/bayroot/slopegraph1.pdf", width=2.5, height=4)
par(mar=c(3,3,1,2))
slopegraph(results$rtt, results$bayroot, names.arg=c('RTT', 'bayroot'), 
           colorize=T, type='l', ylab="Root mean square error", shim=0,
           cex.axis=0.8, mgp=c(2,0.5,0))
dev.off()

