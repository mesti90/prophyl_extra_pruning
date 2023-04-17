# This script takes a phylogenetic tree and attempts to root it using Midpoint
# Rooting.

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
  test_dir <- "~/Methods/prophyl-tests/test-root_tree_mad"
  project_dir <- "~/Methods/prophyl"
  tree_path <- paste0(test_dir, "/sample_tree.tre")
}

# read tree
tree <- ape::read.tree(tree_path)

# if tree is rooted, unroot
if (ape::is.rooted(tree)) {
  tree <- ape::unroot(tree)
}

# root tree
rooted_tree <- phytools::midpoint.root(tree)

# export rooted tree
if (!interactive()) {
  # export rooted tree object
  saveRDS(rooted_tree, "rooted_tree_mp.rds")
  # export newick tree
  writeLines(
    rooted_tree[[1]], 
    con = "rooted_tree_mp.tre"
  )
}

# end logging
if (!interactive()) {
  sink(con)
}
