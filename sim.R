require(twt, quietly=TRUE)
#> 
#> Attaching package: 'ggfree'
#> The following object is masked from 'package:ape':
#> 
#>     unroot
# locate the YAML file that specifies a compartmental SI model
set.seed(1234)
path <- '~/1_projects/bayroot/Simulated_integration.yaml'
model <- Model$new(yaml.load_file(path))

# run an outer tree simulation
outer <- sim.outer.tree(model)
# display the population trajectories
plot(outer, type='s')

inner <- sim.inner.tree(outer)
# display the inner tree annotated with transmission events (points)
plot(inner)

events <- outer$get.eventlog()$get.all.events()
lookup <- list()
lookup[[events$compartment2[nrow(events)]]] <- 1

for (row in nrow(events):4) {
  lookup[[events$compartment1[row]]] <- lookup[[events$compartment2[row]]] + 1
}


for (row in nrow(events):4) {
  events$time[row] <- lookup[[events$compartment2[row]]]
}

string <- "(((6:2,L2:2,4:0):1,8:1,L1:0):1,(3:1,(5:0,L3:1):1):1,(2:0,7:4):2,1:0);"
tr <- read.tree(text=string)
plot(tr)
