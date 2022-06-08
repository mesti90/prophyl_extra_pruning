rm(list=ls())
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
combined_file <- paste0(args[1], "/combined_ancestral_states.tab")
value <- args[2]

combined <- read.csv(combined_file, sep = "\t")
index <- which(names(combined) == value)
names(combined)[index] <- "target"
combined$count <- 1

tidy <- combined %>% complete(node, target, fill = list(count = 0))
tidy <- reshape2::dcast(tidy,node~target, value.var = "count")

tnames <- names(tidy)
tidy <- cbind(
  as.data.frame(tidy[, 1]),
  as.data.frame(t(apply(tidy[,-1], 1, function(x) x/sum(x))))
)
names(tidy) <- tnames

#saveRDS(tidy, paste0(value, "_marginals.rds"))

saveRDS(tidy, "marginals.rds")