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

if (ncol(tidy) > 2) {
  tidy <- cbind(
    as.data.frame(tidy[, 1]),
    as.data.frame(t(apply(tidy[,-1], 1, function(x) x/sum(x))))
  )
} 
names(tidy) <- tnames

tidy <- dplyr::rename(tidy, label = node)

if (ncol(tidy) > 2) {
  ancestral_states <- data.frame(
    label = tidy$label,
    group = value,
    value = apply(tidy[,-1], 1, function(x) {
      paste(names(tidy)[which(x > 0)+1], collapse = "|")
    })
  )
} else {
  ancestral_states <- data.frame(
    label = tidy$label,
    group = value,
    value = names(tidy)[2]
  )
}


write.table(
  ancestral_states,
  file = paste0(value, ".tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#saveRDS(ancestral_states, "ancestral_states.rds")