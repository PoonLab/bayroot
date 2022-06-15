# bayroot: a Bayesian approach to root-to-tip regression

`bayroot` is a simple R package that provides functions for a Bayesian approach to root-to-tip regression.

## Description
Root-to-tip regression (RTT) makes use of the [molecular clock](https://en.wikipedia.org/wiki/Molecular_clock) assumption that the rate of evolution is roughly constant over time.  Starting with a phylogenetic tree relating genetic sequences that have been observed at different points in time, this method regresses the divergence of each sequence from an inferred common ancestor at the root of the tree against the time of sampling.  If the molecular clock assumption holds, then we expect to see a linear increase in divergence over time.  The standard approach to fitting this linear model is to minimize a cost function such as least squares or [R-squared](https://en.wikipedia.org/wiki/Coefficient_of_determination) (Rambaut *et al* 2016).

A more complex approach to this problem is to use Bayesian methods to sample trees from the posterior distribution defined by the sequences, sampling dates, model of evolution, and prior belief.  However, there are applications where a full Bayesian treatment may be difficult to compute, such as for large numbers of sequences, or when many sequences have unknown sampling dates to be estimated.  Therefore we developed `bayroot` to provide a solution that is intermediate between these two extremes.

## Usage
```R
# load example tree (51 influenza HA sequences from 2001-2006) from Russell et al. 2008)
phy <- read.tree("data/h3n2.nwk")

# first show root-to-tip regression using ape::rtt
tip.dates <- get.dates(phy, format="%Y")
rooted <- rtt(phy, as.integer(tip.dates))
div <- node.depth.edgelength(rooted)[1:Ntip(rooted)]
plot(tip.dates, div); abline(lm(div ~ tip.dates))
```
<img src="https://user-images.githubusercontent.com/1109328/173941968-d2b5f464-8369-4746-90b0-a2ebd7b2c427.png" width="300px"/>

```R
# now carry out Bayesian sampling
settings <- list(seq.len=987, format="%Y", 
                 mindate=as.Date("1990-01-01"), maxdate=as.Date("2000-01-01"),
                 meanlog=-5, sdlog=2,
                 root.delta=0.01, date.sd=10, rate.delta=1e-3)
params <- list(phy=rooted, rate=0.01, origin=as.Date("1995-01-01"))  # takes about a minute
res <- bayroot(nstep=1e4, params=params, settings=settings)
plot(res, burnin=100)  # calls a generic S3 method
```
<img src="https://user-images.githubusercontent.com/1109328/173942067-934a648f-c112-480d-b684-135bba19d56c.png" width="500px"/>


## Dependencies

`bayroot` was developed in R version 3.6.2.  It requires the following packages:
* [phytools](https://github.com/liamrevell/phytools) for the `reroot` function
* [msm](https://cran.r-project.org/web/packages/msm/index.html) for the truncated normal distribution function (`rtnorm`)
* [ggfree](https://github.com/ArtPoon/ggfree) for analyzing and plotting trees
* [twt](https://github.com/PoonLab/twt) (*optional*, for simulating test data)


### References
* Rambaut, A., Lam, T. T., Max Carvalho, L., & Pybus, O. G. (2016). Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). Virus evolution, 2(1), vew007.
* Russell, C. A., Jones, T. C., Barr, I. G., Cox, N. J., Garten, R. J., Gregory, V., ... & Smith, D. J. (2008). The global circulation of seasonal influenza A (H3N2) viruses. Science, 320(5874), 340-346.
