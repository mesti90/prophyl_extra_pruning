# This script is used for collecting dated trees from multiple runs into a
# a list of trees. In this project it is used for collecting dated trees that 
# were built for each relevant MLST to a list of trees and using this list in
# downstream analysis.

# directory which contains a number of trees
dirname <- "trees"
jobs <- dir(dirname)

trees <- list()

for (i in seq_along(jobs)) {
  tree <- try(readRDS(paste0(dirname, "/", jobs[i])))
  if (inherits(tree, "try-error")) {
    stop("Input directory must contain rds files and only rds files.")
  }
  testthat::expect_true("treedater" %in% class(tree))
  trees[[i]] <- tree
  names(trees)[i] <- strsplit(jobs[i], "_|\\.")[[1]][3]
}

saveRDS(trees, file = "all_dated_trees.rds")
