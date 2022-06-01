require(ggplot2)
require(chemCal)

# set the file location
tree.file <- '~/1_projects/bayroot/zambia/ZM1165M'

# read tree
tree <- read.tree(paste(tree.file, ".fa.hyp.treefile", sep = ""))
# Numeric date origin
time.since.date <- as.Date('2021-01-01') - as.numeric(as.Date('2021-01-01'))

get.dates <- function(phy, delimiter='_', pos=-1, format='%Y-%m-%d') {
  dt <- sapply(phy$tip.label, function(x) {
    tokens <- strsplit(x, delimiter)[[1]]
    if (pos == -1) {
      return(tokens[length(tokens)])
    }
    else {
      return(tokens[pos])
    }
  })
  as.Date(dt, format=format)
}

data <- data.frame(label=tree$tip.label, 
                   type=gsub(".*((RNA)|(DNA)).*", "\\1", tree$tip.label), 
                   censored=as.numeric(grepl(".*DNA.*", tree$tip.labe)), 
                   date=get.dates(tree),
                   num.date=as.numeric(get.dates(tree)),
                   stringsAsFactors=F)

# remove zero length branches if they exist
if (any(tree$edge.length < 1e-7)) {
  warning("Tiny length branches detected. This could be caused by duplicate sequences.")
}

# root tree
# censor DNA sequences
data$date[data$censored == 1] <- NA

rooted.tree <- rtt(tree, data$num.date, opt.tol = 1e-16)

# linear regression
data$dist <- node.depth.edgelength(rooted.tree)[1:Ntip(rooted.tree)]

model <- lm(dist ~ num.date, data = data, subset = censored == 0)

est.date <- as.data.frame(do.call(
  rbind,
  lapply(data$dist, function(x)
    unlist(inverse.predict(model, x))
  )
))
# Convert prediction to date
est.date$prediction.date <- time.since.date + est.date$Prediction
est.date$CI_low <- time.since.date + est.date$`Confidence Limits1`
est.date$CI_high <- time.since.date + est.date$`Confidence Limits2`
est.date$label <- data$label
est.date$censored <- data$censored

root.prediction <- inverse.predict(model, 0)
root.date <- data.frame('root.date' = time.since.date + root.prediction$Prediction,
                        'CI_low' = time.since.date + root.prediction$`Confidence Limits`[1],
                        'CI_high' = time.since.date + root.prediction$`Confidence Limits`[2])

write.table(est.date, paste(tree.file, ".rtt.csv", sep = ""), row.names = F)

ggplot() + geom_boxplot(aes(y = prediction.date, 
                            min = CI_low,
                            max = CI_high,
                            middle = prediction.date, 
                            lower = CI_low, 
                            upper = CI_high, 
                            x = label), subset(est.date, subset = censored == 1), 
                        stat = "identity") +
  labs(x = "Sequence", y = "Integration date") +
  theme(axis.text.x = element_text(angle = 90))

