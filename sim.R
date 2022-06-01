require(twt, quietly=TRUE)

setwd("~/git/bayroot/")
model <- Model$new(yaml.load_file('latent1.yaml'))

outer <- sim.outer.tree(model)
plot(outer, type='s')  # to check population trends



# colour branch below (towards root) internal node based on label
# US_Active_28__US_Active_23 means transmission from US_Active_28 to US_Active_23
phy <- as.phylo(outer)
L <- tree.layout(phy)
pch <- list('transmission'=17, 'transition'=1, 'tip'=19)

plot(L, type='n') #, label='b', srt=30, cex=0.5)
lines(L, col=ifelse(L$edges$from.type=='Active', 'red', 'blue'))
points(L$edges$x1, L$edges$y1, cex=1,
       pch=as.integer(pch[L$edges$event.type]))


# zero-out branch lengths for latent compartments
set.seed(1)
for (i in 1:10) {
  print(i)
  outer <- sim.outer.tree(model)
  phy <- as.phylo(outer)
  
  # TODO: label tips wrt. sample times
  
  int.times <- get.integration.times(phy)
  
  write.tree(phy, file=paste("data/latent1", 'orig', i, 'nwk', sep='.'))
  phy$edge.length[phy$from.type=='Latent'] <- 0
  write.tree(phy, file=paste("data/latent1", 'cens', i, 'nwk', sep='.'))
}


get.integration.times <- function(phy) {
  L <- tree.layout(phy)  # yields nice data structures
  df <- L$edges[nrow(L$edges):1, ]
  ltips <- which(df$isTip & df$to.type=='Latent')
  result <- rep(NA, length(ltips))
  for (j in 1:length(ltips)) {
    ltip <- ltips[j]
    is.latent <- TRUE
    # traverse down tree towards root
    for (i in (ltip+1):nrow(df)) {
      row <- df[i,]
      if (row$event.type=='transmission') {
        if (is.latent) {
          if (row$from.type=='Active' & row$to.type=='Latent') {
            result[j] <- row$x1
            break
          }
        }
        else {
          # has transitioned to active
          if (row$from.type=='Active' & row$to.type=='Active') {
            result[j] <- row$x1
            break
          }
        }
      } 
      else if (row$event.type=='transition') {
        if (row$from.type=='Active' & row$to.type=='Latent') {
          is.latent <- FALSE
        }
        else if (row$from.type=='Latent' & row$to.type=='Active') {
          is.latent <- TRUE
        }
      }
    }
  }
  names(result) <- L$nodes$label[df$child[ltips]]
  result
}


