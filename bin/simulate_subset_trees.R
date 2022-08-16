rm(list=ls())
library(treedater)

args <- commandArgs(trailingOnly = TRUE)
subset_id <- args[1]
dated_tree <- readRDS(args[2])
nreps <- args[3]
ncpu <- args[4]

simtrees <- treedater::parboot(
  dated_tree,
  nreps = nreps,
  ncpu = ncpu,
  quiet = FALSE
)

for (i in 1:length(simtrees$trees)) {
  simtrees$trees[[i]]$tip.label <- unname(sapply(
    simtrees$trees[[i]]$tip.label, function(x) {
    strsplit(x, "\\|")[[1]][1]
  }))
}

saveRDS(simtrees, file = paste0(subset_id, ".rds"))