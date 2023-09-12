rm(list=ls())

library(treedater)

if(!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  subset_id <- args[1]
  dated_tree_path <- args[2]
  nreps <- as.numeric(args[3])
  ncpu <- as.numeric(args[4])
  launchDir <- args[5]
  resdir <- args[6]
} else {
  test_dir <- "~/Methods/prophyl-tests/test-simulate_subset_trees"
  subset_id <- "test"
  dated_tree_path <- paste0(
    test_dir, "/results/date_subset_tree/subsample_0001/dated_tree.rds")
  nreps <- 1
  ncpu <- 10
  launchDir <- test_dir
  resdir <- "~/Methods/prophyl-tests/test-simulate_subset_trees/results"
}

dated_tree <- readRDS(dated_tree_path)

if (nreps == 1) {
  simtrees <- list(trees = list(dated_tree))
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

if (!interactive()) {
  # export path to the simulated tree as txt
  write.table(
    paste0(resdir, "/simulate_subset_trees/", subset_id, ".rds"),
    file = paste0(subset_id, ".txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  # export simulated tree
  saveRDS(simtrees, file = paste0(subset_id, ".rds"))
}
