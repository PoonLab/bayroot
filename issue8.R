setwd('~/git/bayroot')
source('bayroot.R')

# screened for hypermutation, aligned (MAFFT) and reconstructed tree (IQTREE)
phy <- read.tree('data/ZM1044M.fa.hyp.treefile')

# parse tip dates
tip.dates <- get.dates(phy, pos=4)
tip.dates[grepl("_DNA_", phy$tip.label)] <- NA  # censor DNA tips

phy <- reroot(phy, which.min(tip.dates))
div <- node.depth.edgelength(phy)[1:Ntip(phy)]


origin <- as.Date('2006-01-26')  # seroconversion midpoint
rate <- 0.0144/365


obj.func <- function(rate) {
  -lf(phy, origin, rate)
}
x <- seq(0.1e-5, 5e-5, length.out=30)
plot(x, sapply(x, obj.func))

res <- optim(rate, obj.func, method="Brent", lower=0, upper=1e-4)
optim.rate <- res$par

plot(tip.dates, div)
max.date <- max(get.dates(phy))
segments(x0=origin, x1=max.date, y0=0, y1=rate*(max.date - origin), lty=2)
segments(x0=origin, x1=max.date, y0=0, y1=optim.rate*(max.date - origin))

L <- as.double(optim.rate*(max.date-origin))  # lambda
x <- seq(0, 2, 0.01)
plot(x, pgamma(x, L, lower=F)*gamma(L)/gamma(x), type='l')



qcpois <- function(q, lambda) {
  # solve for outcome x such that CDF = q
  obj.func2 <- function(x) {
    
  }
}
