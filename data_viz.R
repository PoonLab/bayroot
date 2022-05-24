require(dplyr)

tree.file <- '~/1_projects/bayroot/zambia/ZM1165M'
rtt_table <- read.table(paste(tree.file, ".rtt.csv", sep = ""), header = T)
bayroot_table <- read.table(paste(tree.file, ".bayroot.pred.table", sep =""), header = T, stringsAsFactors = F)

rtt_table$label <- gsub("(MT[0-9][0-9][0-9][0-9][0-9][0-9]).*([0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9])", "\\1", rtt_table$label)
rtt_table$min <- rtt_table$CI_low
rtt_table$max <- rtt_table$CI_high
rtt_table$prediction.date <- as.Date(rtt_table$prediction.date, origin = "1970-01-01")
rtt_table$Prediction <- as.Date(rtt_table$Prediction, origin = "1970-01-01")
rtt_table$method <- "rtt"

m_bayroot_table <- reshape2::melt(bayroot_table)
m_bayroot_table$label <- as.character(m_bayroot_table$variable)
m_bayroot_table$value <- as.Date(m_bayroot_table$value, origin = "1970-01-01")
m_bayroot_table$label <- gsub("(MT[0-9][0-9][0-9][0-9][0-9][0-9]).*([0-9][0-9][0-9][0-9].[0-9][0-9].[0-9][0-9])", "\\1", m_bayroot_table$label)

plot_data <- m_bayroot_table %>% group_by(variable) %>%
  dplyr::summarise(Prediction = median(value),
                   Standard.Error = sd(value),
                   Confidence = NA,
                   Confidence.Limits1 = NA,
                   Confidence.Limits2 = NA,
                   prediction.date = quantile(value, probs = 0.5, type = 1),
                   CI_low = quantile(value, probs = 0.25, type = 1),
                   CI_high = quantile(value, probs = 0.75, type = 1),
                   label = label[1],
                   censored = 0,
                   min = min(value),
                   max = max(value),
                   method = "bayroot")

plot_data <- rbind(plot_data[,-1], subset(rtt_table, censored == 1))

ggplot() + geom_boxplot(aes(y = Prediction, 
                            min = min,
                            max = max,
                            middle = prediction.date, 
                            lower = CI_low, 
                            upper = CI_high, 
                            x = label,
                            col = method), 
                        subset(plot_data, label %in% sample(plot_data$label, 3)), stat = "identity") +
  labs(x = "Sequence") +
  coord_flip() +
  theme(axis.text.y = element_blank(), text = element_text(size = 28), legend.position = "none")

trees <- list.files("~/1_projects/bayroot/zambia/", ".fa.hyp.treefile", full.names = T)

divergence <- NULL
pat <- NULL

for (i in trees) 
{
  tree <- read.tree(i)
  rooted.tree <- rtt(tree, as.numeric(get.dates(tree)), opt.tol = 1e-16)
  DNA_sequences <- grepl("(DNA)", rooted.tree$tip.label)

  divergence <- c(divergence, 
                  node.depth.edgelength(rooted.tree)[1:Ntip(rooted.tree)][DNA_sequences])
  pat <- c(pat, 
           rep(gsub("/home/rouxcil/1_projects/bayroot/zambia//(.*).fa.hyp.treefile", "\\1", i),
           sum(DNA_sequences)))
}

divergence_data <- data.frame("Participant" = pat, divergence = divergence)

ggplot() + 
  geom_density(aes(x = divergence, fill = Participant, col = Participant), alpha = 0.08, divergence_data) +
  theme(text = element_text(size = 16)) +
  facet_grid(Participant~., scales = "free")

pdf <- function(t, y, rate, t0, tmax) {
  # probability of integration time (t) given divergence (y)
  L <- as.double(rate*(t-t0))
  rate * L^y * exp(-L) / gammainc(y+1, as.double(rate*(tmax-t0)))
}

sims <- readRDS(paste(tree.file, ".bayroot.RDS", sep = ""))

rate <- median(sims$log$rate)
origin <- median(sims$log$origin) 
tree <- root(read.tree(paste(tree.file, ".fa.hyp.treefile", sep = "")), 1, r = T)

tip.dates <- get.dates(tree)[grepl("RNA", tree$tip.label)]
max.date <- max(tip.dates, na.rm=T)
x <- seq(origin, max.date, length.out=100)

y1 <- pdf(x, 0.04017437, rate, origin, max.date)
y2 <- pdf(x, 0.1495772, rate, origin, max.date)
y3 <- pdf(x, 0.2169257, rate, origin, max.date)

h <- 1/as.integer(max.date - origin)
par(mfrow=c(1,3), cex.lab=1.2)
plot(x, y1, type='l', xlab='Sampling date', ylab='Probability density',
     main=0.04)
abline(h=h, lty=2)
plot(x, y2, type='l', col='red', xlab='Sampling date', 
     ylab='Probability density', main=0.15)
abline(h=h, lty=2)
plot(x, y3, type='l', col='blue', xlab='Sampling date', 
     ylab='Probability density', main=0.22)
abline(h=h, lty=2)

