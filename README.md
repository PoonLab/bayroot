# bayroot: a Bayesian approach to root-to-tip regression

`bayroot` is a simple R package that provides functions for a Bayesian approach to root-to-tip regression.

## Description
Root-to-tip regression (RTT) makes use of the [molecular clock](https://en.wikipedia.org/wiki/Molecular_clock) assumption that the rate of evolution is roughly constant over time.  Starting with a phylogenetic tree relating genetic sequences that have been observed at different points in time, this method regresses the divergence of each sequence from an inferred common ancestor at the root of the tree against the time of sampling.  If the molecular clock assumption holds, then we expect to see a linear increase in divergence over time.  The standard approach to fitting this linear model is to minimize a cost function such as least squares or [R-squared](https://en.wikipedia.org/wiki/Coefficient_of_determination) (Rambaut *et al* 2016).

A more complex approach to this problem is to use Bayesian methods to sample trees from the posterior distribution defined by the sequences, sampling dates, model of evolution, and prior belief.  However, there are applications where a full Bayesian treatment may be difficult to compute, such as for large numbers of sequences, or when many sequences have unknown sampling dates to be estimated.  Therefore we developed `bayroot` to provide a solution that is intermediate between these two extremes.

## Usage

### Fitting a root-to-tip regression
```R
# load example tree (101 influenza HA sequences from 2001-2006) from Russell et al. 2008)
phy <- read.tree("data/h3n2.nwk")

# first show root-to-tip regression using ape::rtt
tip.dates <- get.dates(phy, format="%d-%b-%Y")
rooted <- rtt(phy, as.integer(tip.dates))
div <- node.depth.edgelength(rooted)[1:Ntip(rooted)]
plot(tip.dates, div, xlim=c(as.Date("1999-01-01"), max(tip.dates)), 
     ylim=c(0, max(div)))
abline(lm(div ~ tip.dates))
```
<img src="https://user-images.githubusercontent.com/1109328/174401842-535dda29-2c62-42db-9cf0-de556dd9acd3.png" width="400px"/>


```R
# now carry out Bayesian sampling, using the RTT tree to initialize the chain
settings <- list(seq.len=987, format="%d-%b-%Y", 
                 mindate=as.Date("1995-01-01"), maxdate=as.Date("2000-01-01"),
                 meanlog=-5, sdlog=2,
                 root.delta=0.03, date.sd=60, rate.delta=0.002)
params <- list(phy=rooted, rate=0.01, origin=as.Date("1995-01-01"))
res <- bayroot(nstep=1e5, skip=100, params=params, settings=settings, echo=T)  # about 8 minutes
plot(res, burnin=100)  # calls a generic S3 method
```
<img src="https://user-images.githubusercontent.com/1109328/174401648-b6691a29-65fc-4104-9b37-a7e75dd1bf04.png" width="600px"/>


```R
# since we can't make a trace log of trees, bayroot provides a plotting option for viewing a specific step in the chain sample
plot(res, settings, step=1000)
```
<img src="https://user-images.githubusercontent.com/1109328/174401959-6c7d857e-6189-4574-bb64-c6887f6978ba.png" width="600px"/>

```R
# another way of looking at the MCMC results
plot(tip.dates, div, xlim=c(as.Date("1999-01-01"), max(tip.dates)), 
      ylim=c(0, max(div)), type='n')
abline(fit)
for (i in seq(100, nrow(res$log), length.out=100)) {
  row <- res$log[i,]
  clock <- row$rate / settings$seq.len
  segments(
    x0=row$origin, x1=max(tip.dates), 
    y0=0, y1=clock*as.integer(max(tip.dates)-row$origin),
    col=rgb(0,0,1,0.2)
  )
}
points(tip.dates, div, pch=21, cex=1, bg='white')
```
<img src="https://user-images.githubusercontent.com/1109328/174402219-072167a9-ee2d-4baa-938e-8214abf6db6d.png" width="400px"/>

### Estimating unknown tip dates

```R
require(ape)

# 50 trees generated from a simulation of cell dynamics
latent1 <- read.tree("data/testdata-latent1.nwk")


```

## Dependencies

`bayroot` was developed in R version 3.6.2.  It requires the following packages:
* [phytools](https://github.com/liamrevell/phytools) for the `reroot` function
* [msm](https://cran.r-project.org/web/packages/msm/index.html) for the truncated normal distribution function (`rtnorm`)
* [ggfree](https://github.com/ArtPoon/ggfree) for analyzing and plotting trees
* [twt](https://github.com/PoonLab/twt) (*optional*, for simulating test data)


### References
* Rambaut, A., Lam, T. T., Max Carvalho, L., & Pybus, O. G. (2016). Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). Virus evolution, 2(1), vew007.
* Russell, C. A., Jones, T. C., Barr, I. G., Cox, N. J., Garten, R. J., Gregory, V., ... & Smith, D. J. (2008). The global circulation of seasonal influenza A (H3N2) viruses. Science, 320(5874), 340-346.
