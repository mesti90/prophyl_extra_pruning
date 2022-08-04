rm(list=ls())
library(treedater)

args <- commandArgs(trailingOnly = TRUE)
dated_tree <- readRDS(args[1])
nreps <- args[2]
ncpu <- args[3]

simulated_trees <- treedater::parboot(
  dated_tree,
  nreps = nreps,
  ncpu = ncpu,
  quiet = FALSE
)

saveRDS(simulated_trees, file = "simulated_trees.rds")