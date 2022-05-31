require(twt, quietly=TRUE)

# locate the YAML file that specifies a compartmental SI model
set.seed(1)

#path <- '~/1_projects/bayroot/Simulated_integration.yaml'
path <- '~/git/bayroot/latent.yaml'
model <- Model$new(yaml.load_file(path))
outer <- sim.outer.tree(model)


plot(outer, type='s')


phy <- as.phylo(outer)
par(xpd=NA)
plot(phy, show.node.label=T, cex=0.5, srt=30)
par(xpd=F)


if (FALSE) {
  par(mar=c(5,5,1,1))
  plot(NA, xlim=c(-20, 0), ylim=c(0, 1e3), xlab="Time", ylab="Count")
  abline(v=-10)
  for (i in 1:10) {
    outer <- sim.outer.tree(model)
    counts <- outer$get.counts()
    y.max <- max(counts$I.Active)
    lines(-counts$time, counts$S.Active, col='blue')
    lines(-counts$time, counts$I.Active, col='red')
    lines(-counts$time, counts$S.Latent, col='dodgerblue')
    lines(-counts$time, counts$I.Latent, col='salmon')
  }
  legend(x=-8, y=900, legend=c("S.Active", "I.Active", "S.Latent", "I.Latent"),
         col=c('blue', 'red', 'dodgerblue', 'salmon'), lwd=2, bty='n')
}

# colour branch below (towards root) internal node based on label
# US_Active_28__US_Active_23 means transmission from US_Active_28 to US_Active_23

L <- tree.layout(phy)

pch <- list('transmission'=17, 'transition'=1, 'tip'=19)

plot(L, type='n') #, label='b', srt=30, cex=0.5)
lines(L, col=ifelse(L$edges$from.type=='Active', 'red', 'blue'))
points(L$edges$x1, L$edges$y1, cex=1,
       pch=as.integer(pch[L$edges$event.type]))


# zero-out branch lengths for latent compartments
require(twt)
setwd('~/git/bayroot')
path <- 'latent1.yaml'
model <- Model$new(yaml.load_file(path))
set.seed(1)
for (i in 1:10) {
  print(i)
  outer <- sim.outer.tree(model)
  phy <- as.phylo(outer)
  write.tree(phy, file=paste("data/latent1", 'orig', i, 'nwk', sep='.'))
  phy$edge.length[phy$from.type=='Latent'] <- 0
  write.tree(phy, file=paste("data/latent1", 'cens', i, 'nwk', sep='.'))
}

#plot(tree.layout(phy))

eventlog <- outer$get.eventlog()
events <- eventlog$get.all.events()

