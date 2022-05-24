source('~/1_projects/bayroot/bayroot.R')

# set the file location
tree.file <- '~/1_projects/bayroot/zambia/ZM1044M'
# read tree
tree <- root(read.tree(paste(tree.file, ".fa.hyp.treefile", sep = "")), 1, r = T)

data <- data.frame(label=tree$tip.label, 
                   type=gsub(".*((RNA)|(DNA)).*", "\\1", tree$tip.label), 
                   censored=as.numeric(grepl(".*DNA.*", tree$tip.labe)), 
                   date=get.dates(tree),
                   num.date=as.numeric(get.dates(tree)),
                   stringsAsFactors=F)

params <- list('phy' = tree, 'origin' = as.Date('1970-01-01'), 
               'rate' = 1.228e-05)
settings <- list('mindate' = as.Date('1960-01-01'), 
                 'maxdate' = as.Date('2020-01-01'), 
                 'meanlog' = 0, 
                 'sdlog' = 1, 
                 'root.delta' = NA,
                 'date.sd' = 1825,
                 'rate.delta' = 1e-04)

test <- mh(10000, params, settings)

saveRDS(test, paste(tree.file, ".bayroot.RDS", sep =""))

test2 <- predict.bayroot(test, data$label[data$censored==1])
write.table(test2, paste(tree.file, ".bayroot.pred.table", sep =""), row.names = F)

bayroot.predict <- reshape2::melt(test2)
bayroot.predict$value <- as.Date(bayroot.predict$value, origin = "1970-01-01")

dodge <- position_dodge(width=-1) 
ggplot() + geom_boxplot(aes(y = prediction.date, 
                            min = CI_low,
                            max = CI_high,
                            middle = prediction.date, 
                            lower = CI_low, 
                            upper = CI_high, 
                            x = label), subset(est.date, subset = censored == 1), stat = "identity") +
  geom_boxplot(aes(y = value, x = variable), data = bayroot.predict, col = "red", position = dodge) +
  theme(axis.text.x = element_text(angle = 90))
  

       