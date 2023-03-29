# This script takes a phylogenetic tree and attempts to root it using Minimum
# Ancestor Deviation (MAD).

rm(list = ls())

# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

# define input
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  project_dir <- args[1]
  tree_path <- args[2]
} else {
  test_dir <- "~/Methods/prophyl-tests/test_root_tree_mad"
  project_dir <- "~/Methods/prophyl"
  tree_path <- paste0(test_dir, "/treeshrink.tre")
}

library(devtools)
load_all(project_dir)

# read tree
tree <- ape::read.tree(tree_path)

# root tree
rooted_tree <- root_mad(tree, output_mode = "full", verbose = TRUE)

# export rooted tree
if (!interactive()) {
  ape::write.tree(
    rooted_tree, 
    file = "rooted_tree.tre"
  )
}

# end logging
if (!interactive()) {
  sink(con)
}
