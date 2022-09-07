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

First, we'll load the *bayroot* functions and a couple of packages, and then import the simulated data:
```R
source("bayroot.R")
require(ape)
require(lubridate)

# 50 trees generated from a simulation of cell dynamics
latent1 <- read.tree("data/testdata-latent1.nwk")
int.times <- read.csv("data/testdata-latent1.csv")
```

Next, we'll extract the first tree and the actual integration dates, and configure our analysis:
```R
phy <- latent1[[1]]  # extract the first replicate
int.times <- int.times[int.times$rep=='1',]

settings <- list(
  seq.len=1233,  # AY772699
  format="%Y-%m-%d",
  # we arbitrarily assigned 2000-01-01 as the actual origin date
  mindate=as.Date("1999-12-01"),  # one month before origin
  maxdate=as.Date("2000-04-01"),  # first sample is 3 months after origin
  meanlog=-5, sdlog=2,  # hyperparameters for lognormal prior on rate
  root.delta=0.01,  # proposal on root location
  date.sd=10,  # days, origin proposal
  rate.delta=0.01,  # proposal on clock rate
  censored = phy$tip.label[grepl("^Latent", phy$tip.label)]
)
```

Next, we're going to convert the sampling times in the tip labels into dates, and identify which tip dates need to be imputed:
```R
tip.dates <- sapply(phy$tip.label, function(x) as.integer(strsplit(x, "_")[[1]][3]))
tip.dates <- as.Date("2000-01-01") + months(tip.dates)
phy$tip.label <- paste(phy$tip.label, tip.dates, sep="_")

# if tip was sampled from latent compartment, censor sample date
censored <- as.Date(ifelse(grepl("Active", phy$tip.label), tip.dates, NA), 
                      origin='1970-01-01')
```

Finally, we'll initialize the model parameters, and then run a chain sample for 20,000 steps, only keeping every 20th step.  I've set `echo=TRUE` to display these samples (this will take a couple of minutes):
```R
# initial model parameters
params <- list(
  phy=reroot(phy, which.min(censored)),  # initially root on earliest tip
  rate=0.1,  # initial clock rate
  origin=min(censored, na.rm=T)-1  # initial date at the root
  )
set.seed(1)
chain <- bayroot(nstep=2e4, skip=20, params=params, settings=settings, echo=T)
```

If we call the generic `plot` method on the object returned by `bayroot()`, R will display a composite set of plots summarizing the posterior traces of our chain sample.  We'll discard the first 2,000 steps as burnin:
```R
plot(chain, burnin=100)
```
<img src="https://user-images.githubusercontent.com/1109328/188780713-72137dc0-7c98-4e63-ab64-720191234ac1.png" width="600px"/>

Now we use this posterior sample to simulate integration dates for each of the censored tips.  Note that for this simulation, we assumed that ART was initiated exactly 11 months post-infection, so we have to pass this information to the `max.date` argument:
```R
pred.dates <- predict(chain, settings, max.date=as.Date("2000-11-01), burnin=100, thin=200)
```

This function returns a list of data frames for each censored tip:
```R
> summary(pred.dates)
                            Length Class      Mode
Latentcomp_5_20_2001-09-01  2      data.frame list
Latentcomp_7_20_2001-09-01  2      data.frame list
Latentcomp_3_20_2001-09-01  2      data.frame list
Latentcomp_10_20_2001-09-01 2      data.frame list
Latentcomp_4_20_2001-09-01  2      data.frame list
Latentcomp_2_20_2001-09-01  2      data.frame list
Latentcomp_8_20_2001-09-01  2      data.frame list
Latentcomp_6_20_2001-09-01  2      data.frame list
Latentcomp_1_20_2001-09-01  2      data.frame list
Latentcomp_9_20_2001-09-01  2      data.frame list
> summary(pred.dates[[1]])
    int.date          div       
 Min.   :11056   Min.   :36.94  
 1st Qu.:11088   1st Qu.:40.13  
 Median :11104   Median :41.29  
 Mean   :11106   Mean   :41.10  
 3rd Qu.:11121   3rd Qu.:42.19  
 Max.   :11169   Max.   :44.45
```

Let's compare these estimates to standard root-to-tip regression using a couple of utility functions that we've provided in `scripts/validate.R`:
```R
source("validate.R")
rt <- root2tip(phy)
get.true.values.1(rt, int.times)
est <- get.estimates.1(pred.dates, rt)
```

and generate a plot summarizing this comparison:
```R
par(mar=c(5,5,1,1))
plot(rt, true.vals=true.vals, xlab="Collection date (months since origin)",
     ylab="Divergence", xlim=c(0, 20), ylim=c(0, max(rt$div)),
     las=1, cex.axis=0.8)
points(est$est, est$div, pch=19, col=rgb(0,0,1,0.5), cex=0.8)
segments(x0=est$lo95, x1=est$hi95, y0=est$div, col=rgb(0,0,1,0.3), lwd=5)
abline(v=10, lty=2)
legend(x=1, y=0.035, legend=c("RTT", "bayroot", "true date"), cex=0.8,
       col=c('red', 'blue', 'black'), pch=c(19, 19, 3), pt.lwd=2, bty='n')
```
<img src="https://user-images.githubusercontent.com/1109328/188785869-ffed8d02-54d8-40a3-90f5-85a61faebe92.png" width="500px"/>
Note the open circles represent RNA sequences used to calibrate the molecular clock.


### References
* Rambaut, A., Lam, T. T., Max Carvalho, L., & Pybus, O. G. (2016). Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). Virus evolution, 2(1), vew007.
* Russell, C. A., Jones, T. C., Barr, I. G., Cox, N. J., Garten, R. J., Gregory, V., ... & Smith, D. J. (2008). The global circulation of seasonal influenza A (H3N2) viruses. Science, 320(5874), 340-346.
