rm(list=ls())
library(treedater)

args <- commandArgs(trailingOnly = TRUE)
subset_id <- args[1]
dated_tree <- readRDS(args[2])
nreps <- as.numeric(args[3])
ncpu <- as.numeric(args[4])
launchDir <- args[5]

if (nreps == 1) {
  simtrees <- dated_tree
} else {
  simtrees <- treedater::parboot(
    dated_tree,
    nreps = nreps,
    ncpu = ncpu,
    quiet = FALSE
  )
}

for (i in 1:length(simtrees$trees)) {
  simtrees$trees[[i]]$tip.label <- unname(sapply(
    simtrees$trees[[i]]$tip.label, function(x) {
    strsplit(x, "\\|")[[1]][1]
  }))
}

write.table(
  paste0(launchDir, "/results/simulate_subset_trees/", subset_id, ".rds"),
  file = paste0(subset_id, ".txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

saveRDS(simtrees, file = paste0(subset_id, ".rds"))