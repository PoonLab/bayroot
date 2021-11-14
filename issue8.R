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

# find optimal rate given fixed origin and root
obj.func <- function(rate) {
  -lf(phy, origin, rate)
}
x <- seq(0.1e-5, 5e-5, length.out=30)
plot(x, sapply(x, obj.func))  # view likelihood function profile

res <- optim(rate, obj.func, method="Brent", lower=0, upper=1e-4)
optim.rate <- res$par

par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(tip.dates, div)
max.date <- max(get.dates(phy))
segments(x0=origin, x1=max.date, y0=0, y1=rate*(max.date - origin), lty=2)
segments(x0=origin, x1=max.date, y0=0, y1=optim.rate*(max.date - origin))


# look at probability distribution at last time point
lambda <- as.double(optim.rate * (max.date-origin))
y <- seq(0, 1, 0.01)  # divergence
barplot(lambda^y * exp(-lambda) / gamma(y+1), names.arg=y, las=2)


# try Equation (5) from Ilienko by numerical integration
require(cubature)
dcpois <- function(x, lambda) {
  f <- function(st) { exp(-sum(st)) * prod(st)^(x-1) * log(st[1]/st[2]) }
  res <- adaptIntegrate(f, lowerLimit=c(0, lambda), upperLimit=c(lambda, Inf))
  res$integral / gamma(x)^2
}

dcpois(2,1)



# Try out cumulative distribution function of continuous Poisson
# According to Ilienko, Annales Univ. Sci. Budapest., Sect. Comp. 39 (2013) 137-147
#    F(x) = \Gamma(x, a) / \Gamma(x)  
# where a is Poisson rate parameter and 
#    \Gamma(x, a) = \int_{a}^{\infty} \exp{-t} t^{x-1} dt
#
# Note this notation is switched w.r.t. R help page for pgamma(x, a):
#    P(a, x) = \int_{x}^{\infty} \exp{-t} t^{a-1} dt
L <- as.double(optim.rate*(max.date-origin))  # lambda
x <- seq(0, 0.1, length.out=100)
plot(x, pgamma(q=L, shape=x, lower=F), type='l')


p <- seq(0,1,length.out=50)
lines(p, qgamma(p=p, shape=L), type='l')

qgamma(p=0.5, shape=L)
points(y=0.5, x=qgamma(p=0.5, shape=L))


############################################################
# Try adjusting divergence by alignment (sequence) length

phy <- read.tree('data/ZM1044M.fa.hyp.treefile')
tip.dates <- get.dates(phy, pos=4)
tip.dates[grepl("_DNA_", phy$tip.label)] <- NA  # censor DNA tips
phy <- reroot(phy, which.min(tip.dates))

# adjust all branch lengths
seqlen <- 9200
phy$edge.length <- phy$edge.length * seqlen

div <- node.depth.edgelength(phy)[1:Ntip(phy)]

origin <- as.Date('2006-01-26')  # seroconversion midpoint
rate <- 0.0144/365 * seqlen

# find optimal rate given fixed origin and root
obj.func <- function(rate) {
  -lf(phy, origin, rate)
}
x <- seq(0.1e-5, 0.5, length.out=100)
plot(x, sapply(x, obj.func))  # view likelihood function profile

res <- optim(rate, obj.func, method="Brent", lower=0, upper=1)
optim.rate <- res$par

par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(tip.dates, div)
max.date <- max(get.dates(phy))
segments(x0=origin, x1=max.date, y0=0, y1=rate*(max.date - origin), lty=2)
segments(x0=origin, x1=max.date, y0=0, y1=optim.rate*(max.date - origin))


# look at probability distribution at last time point
lambda <- as.double(optim.rate * (max.date-origin))
y <- seq(0, 1e3, length.out=100)  # divergence
barplot(exp(y*log(lambda) -lambda - lgamma(y+1)), names.arg=round(y), las=2,
        cex.names=0.8)


