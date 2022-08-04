rm(list=ls())
library(treedater)

args <- commandArgs(trailingOnly = TRUE)
dated_tree <- readRDS(args[1])
nreps <- args[2]
ncpu <- args[3]

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

saveRDS(simtrees, file = "simtrees.rds")