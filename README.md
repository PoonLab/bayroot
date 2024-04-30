# bayroot: a Bayesian approach to root-to-tip regression

`bayroot` is a simple R package that provides functions for a Bayesian approach to root-to-tip regression.

## Description
Root-to-tip regression (RTT) makes use of the [molecular clock](https://en.wikipedia.org/wiki/Molecular_clock) assumption that the rate of evolution is roughly constant over time.  Starting with a phylogenetic tree relating genetic sequences that have been observed at different points in time, this method regresses the divergence of each sequence from an inferred common ancestor at the root of the tree against the time of sampling.  If the molecular clock assumption holds, then we expect to see a linear increase in divergence over time.  The standard approach to fitting this linear model is to minimize a cost function such as least squares or [R-squared](https://en.wikipedia.org/wiki/Coefficient_of_determination) (Rambaut *et al* 2016).

A more complex approach to this problem is to use Bayesian methods to sample trees from the posterior distribution defined by the sequences, sampling dates, model of evolution, and prior belief.  However, there are applications where a full Bayesian treatment may be difficult to compute, such as for large numbers of sequences, or when many sequences have unknown sampling dates to be estimated.  Therefore we developed `bayroot` to provide a solution that is intermediate between these two extremes.


## Dependencies

`bayroot` was developed in R version 3.6.2.  It requires the following packages:
* [phytools](https://github.com/liamrevell/phytools) for the `reroot` function
* [msm](https://cran.r-project.org/web/packages/msm/index.html) for the truncated normal distribution function (`rtnorm`)
* [ggfree](https://github.com/ArtPoon/ggfree) for analyzing and plotting trees
* [chemCal](https://cran.r-project.org/web/packages/chemCal/index.html) for `inverse.predict` function
* [twt](https://github.com/PoonLab/twt) (*optional*, for simulating test data)


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

*bayroot* was actually developed for our study of the latent HIV reservoir.  HIV converts its RNA genome into DNA and then integrates that DNA into the host genome.  In a small number of host cells, that HIV DNA can remain dormant for years until the cell becomes "reactivated" and starts to produce new virus copies.  We are interested in reconstructing *when* copies of HIV in this latent reservoir became integrated.

Since the HIV DNA basically stops evolving upon integration, this problem is equivalent to imputing missing dates on some tips of a time-scaled tree.  A Bayesian approach is well-suited to this problem!  We implemented a compartmental model in the R package [twt](https://github.com/PoonLab/twt) to simulate trees relating a mixture of actively evolving HIV RNA and integrated HIV DNA.  Let's reconstruct the integration dates for one of these simulations.

First, we'll load the *bayroot* functions and a couple of packages, and then import the simulated data.
We wrote some utility functions that were designed to operate on files containing only one tree and integration times for one replicate, but we didn't want to upload hundreds of these files to the repository.  Instead we're going to extract one replicate and generate the expected file formats:
```R
source("bayroot.R")
require(ape)
require(lubridate)

# 50 trees generated from a simulation of cell dynamics
trees <- read.tree("data/testdata-latent2.nwk")
times <- read.csv("data/testdata-latent2.csv")

# export one of the replicates to new files
rep <- 1
tf <- "data/example.nwk"
write.tree(trees[[rep]], file=tf)

cf <- "data/example.csv"
int.times <- times$int.time[times$rep==rep]
names(int.times) <- times$sample[times$rep==rep]
write.csv(as.data.frame(int.times), file=cf)
```

Next, we'll configure our analysis:
```R
phy <- read.tree(tf)
settings <- list(
  seq.len=1233,  # AY772699
  format="%Y-%m-%d",
  # we arbitrarily assigned 2000-01-01 as the actual origin date
  mindate=as.Date("1999-12-01"),  # one month before origin
  maxdate=as.Date("2000-02-01"),  # one month after origin
  meanlog=-5, sdlog=2,  # hyperparameters for lognormal prior on rate
  root.delta=0.01,  # proposal on root location
  date.sd=10,  # days, origin proposal
  rate.delta=0.01,  # proposal on clock rate
  censored = phy$tip.label[grepl("^Latent", phy$tip.label)]
)
```

This utility function will initialize the model parameters and then run a chain sample for 20,000 steps, only keeping every 20th step.
If you want to see the chain being updated in real time, you can add an argument `echo=TRUE` - note this will generate a lot of console output!
Otherwise just wait a couple of minutes for the chain sample to complete.
```R
source("scripts/validate.R")  # loads package chemCal
set.seed(1)  # make this reproducible
res <- fit.bayroot(tf, cf, settings=settings, nstep=2e4, skip=20, 
                   max.date=as.Date("2001-04-01"))
```
Note we set the maximum integration date to the above date because ART was initiated 15 months post-infection in this set of simulations.

Here I'm displaying the contents of the list returned from `fit.bayroot`:
```
> summary(res)
           Length Class   Mode   
phy         5     phylo   list   
pred.dates 10     -none-  list   
tip.dates  40     Date    numeric
censored   40     Date    numeric
chain       2     bayroot list   
```
* `phy` is just the original tree
* `pred.dates` is a list containing data frames of integration date estimates for every censored tip
* `tip.dates` is a vector of sampling times associated with each tip as Date objects, including censored tips
* `censored` replaces Date values in `tip.dates` with missing values (`NA`) for censored tips
* `chain` is a custom [S3](http://adv-r.had.co.nz/S3.html) object of class `bayroot` that contains our chain samples of model parameters and rooted trees from the posterior distribution

If we call the generic `plot` method on the `bayroot` object, R will display a composite set of plots summarizing the posterior traces of our chain sample.  We'll discard the first 2,000 steps as burnin:
```R
plot(res$chain, burnin=100)
```
<img src="https://user-images.githubusercontent.com/1109328/189016258-186c1117-5821-450e-ab12-e0f9cf98cf19.png" width="600px"/>


Let's compare our posterior sample of integration dates to standard root-to-tip regression using a couple of utility functions that we've provided in `scripts/validate.R`:
```R
rt <- root2tip(phy)
true.vals <- get.true.values(rt, cf)
est <- get.estimates(res, rt)
```
and generate a plot summarizing this comparison:
```R
par(mar=c(5,5,1,1))
plot(rt, true.vals=true.vals, xlab="Time since infection",
     ylab="Divergence", xlim=c(0, 20), ylim=c(0, max(rt$div)),
     las=1, cex.axis=0.8)
points(est$est, est$div, pch=19, col=rgb(0,0,1,0.5), cex=0.8)
segments(x0=est$lo95, x1=est$hi95, y0=est$div, col=rgb(0,0,1,0.3), lwd=5)
abline(v=15, lty=2)
legend(x=0, y=0.035, legend=c("RTT", "bayroot", "true date"), cex=0.8,
       col=c('red', 'blue', 'black'), pch=c(19, 19, 3), pt.lwd=2, bty='n')
```
<img src="https://user-images.githubusercontent.com/1109328/189016331-01a1bfa8-c086-4b86-a455-60f1ad13aa11.png" width="500px"/>
Note the open circles represent RNA sequences used to calibrate the molecular clock.


### References
* Rambaut, A., Lam, T. T., Max Carvalho, L., & Pybus, O. G. (2016). Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). Virus evolution, 2(1), vew007.
* Russell, C. A., Jones, T. C., Barr, I. G., Cox, N. J., Garten, R. J., Gregory, V., ... & Smith, D. J. (2008). The global circulation of seasonal influenza A (H3N2) viruses. Science, 320(5874), 340-346.
