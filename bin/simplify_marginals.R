rm(list=ls())
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)

combined <- read.csv(args[1], sep = "\t")
combined$count <- 1

tidy <- combined %>% complete(node, country, fill = list(count = 0))
tidy <- reshape2::dcast(tidy,node~country, value.var = "count")

tnames <- names(tidy)
tidy <- cbind(
  as.data.frame(tidy[, 1]),
  as.data.frame(t(apply(tidy[,-1], 1, function(x) x/sum(x))))
)
names(tidy) <- tnames

saveRDS(tidy, "simplified_marginals.rds")
