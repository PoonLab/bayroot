# work through test case ZM1044M - https://doi.org/10.1371/journal.ppat.1008378

setwd('~/git/bayroot')
source('bayroot.R')

# screened for hypermutation, aligned (MAFFT) and reconstructed tree (IQTREE)
phy <- read.tree('data/ZM1044M.fa.hyp.treefile')
#phy <- midpoint.root(phy)

# parse tip dates
tip.dates <- get.dates(phy, pos=4)
tip.dates[grepl("_DNA_", phy$tip.label)] <- NA  # censor DNA tips

phy <- reroot(phy, which.min(tip.dates))


pdfunc <- function(t, y, rate, t0, tmax) {
  # probability of integration time (t) given divergence (y)
  # note uniform prior on t cancels out
  L <- as.double(rate*(t-t0))
  rate * L^y * exp(-L) / inc.gamma(y+1, as.double(rate*(tmax-t0)))
}

origin <- as.Date("2006-01-26")  # sero midpoint
min.date <- min(tip.dates, na.rm=T)
max.date <- max(tip.dates, na.rm=T)

x <- seq(origin, max.date, length.out=100)

y1 <- pdfunc(x, 0.001, 1e-5, min.date, max.date)
y2 <- pdfunc(x, 0.01, 1e-5, min.date, max.date)
y3 <- pdfunc(x, 0.1, 1e-5, min.date, max.date)

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




