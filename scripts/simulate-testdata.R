# This script was used to simulate trees from a compartmental model
# of the within-host dynamics of actively- and latently-infected cells

require(twt, quietly=TRUE)

setwd("~/git/bayroot/")
prefix <- "latent2"
model <- Model$new(yaml.load_file(paste0(prefix, '.yaml')))

if (FALSE) {
  set.seed(1)
  outer <- sim.outer.tree(model)
  plot(outer, type='s')  # to check population trends
  
  # colour branch below (towards root) internal node based on label
  # US_Active_28__US_Active_23 means transmission from US_Active_28 to US_Active_23
  phy <- as.phylo(outer)
  L <- tree.layout(phy)
  pch <- list('transmission'=17, 'transition'=1, 'tip'=19)
  
  #pdf("issue15.pdf", width=6, height=6)
  plot(L, type='n') #, label='b', srt=30, cex=0.5)
  lines(L, col=ifelse(L$edges$from.type=='Active', 'red', 'blue'))
  points(L$edges$x1, L$edges$y1, cex=1,
         pch=as.integer(pch[L$edges$event.type]))
  #dev.off()  
}


get.integration.times <- function(phy) {
  L <- tree.layout(phy)  # yields nice data structures
  df <- L$edges[nrow(L$edges):1, ]
  max.depth <- max(df$x1)
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
            result[j] <- max.depth-row$x1
            break
          }
        }
        else {
          # has transitioned to active
          if (row$from.type=='Active' & row$to.type=='Active') {
            result[j] <- max.depth-row$x1
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


# get sampling times (in forward time)
fixed.sampl <- model$get.initial.conds()$originTime - model$get.fixed.samplings()
names(fixed.sampl) <- gsub("__.+$", "", names(fixed.sampl))

set.seed(1)
for (i in 1:50) {
  print(i)
  outer <- sim.outer.tree(model)
  phy <- as.phylo(outer)
  
  int.times <- get.integration.times(phy)
  write.csv(as.data.frame(int.times), 
            file=paste0("data/", prefix, ".", i, '.times.csv'))
  
  # append sample times (forward) to tip names
  idx <- match(phy$tip.label, names(fixed.sampl))
  phy$tip.label <- paste(phy$tip.label, fixed.sampl[idx], sep="_")
  write.tree(phy, file=paste0("data/", prefix, ".", i, '.orig.nwk'))
  
  # zero-out branch lengths for latent compartments, for simulating evolution
  phy$edge.length[phy$from.type=='Latent'] <- 0
  write.tree(phy, file=paste0("data/", prefix, ".", i, '.cens.nwk'))
}


# generate a collated CSV
files <- Sys.glob("data/latent2.*.times.csv")
int.times <- {}
for (f in files) {
  # ordering should be the same as testdata-latent1.nwk
  temp <- read.csv(f)
  index <- gsub("data/latent[12]\\.([0-9]+)\\.times\\.csv", "\\1", f)
  temp$rep <- index
  names(temp) <- c("sample", "int.time", "rep")
  int.times <- rbind(int.times, temp)
}
write.csv(int.times, file="data/latent2.int-times.csv", row.names=F, quote=F)
