library(dplyr)
library(optparse)
library(TreeTools)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("-l", "--launch_dir"),
    type = "character",
    help = "Directory from which the pipeline was launched."
  ),
  make_option(
    c("-t", "--tree"),
    type = "character",
    help = "A dated tree in rds format."
  ),
  make_option(
    c("-d", "--duplicates"),
    type = "character",
    help = "A text file which contains tip labels for identical tips."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/prophyl",
    launch_dir = getwd(),
    tree = "final_dated_tree.rds",
    duplicates = "duplicates.txt"
  )
}

library(devtools)
load_all(args$project_dir)

# set cache directory for R.cache
R.cache::setCacheRootPath(paste0(
  args$launch_dir,
  ".cache/R/R.cache"
))

# import tree
tree <- readRDS(args$tree)

# import duplicates
if (grepl("\\.txt$", args$duplicates)) {
  duplicates <- parse_duplicates(args$duplicates)
} else if (grepl("\\.rds", args$duplicates)) {
  duplicates <- readRDS(args$duplicates)
}

if (length(duplicates) > 0) {
  
  # add duplicates to tree
  for (i in seq_along(duplicates)) {
    for (j in 1:length(duplicates[[i]])) {
      if (duplicates[[i]][j] != names(duplicates)[i]) {
        tree <- TreeTools::AddTip(
          tree,
          where = names(duplicates)[i],
          label = duplicates[[i]][j],
          edgeLength = 0
        )
      }
    }
  }
  
  # remove tips which were used as reference for duplicates but were not in tree
  # this is relevant when tree is a subset tree.
  for (i in seq_along(duplicates)) {
    # if the name of the list entry (reference) is not identical with its first
    # element then the reference was not part of the subset and should be
    # removed. 
    if (names(duplicates)[i] != duplicates[[i]][1]) {
      tree <- ape::drop.tip(tree, names(duplicates)[i])
    }
  }
  
  # consistency check
  duplicate_tips <- duplicates %>% unlist() %>% unique()
  testthat::expect_true(all(duplicate_tips %in% tree$tip.label))
  
}

# Rename internal nodes. This is needed because TreeTools::AddTips() does not
# alter node.label when adding tips. Raised an issue here:
# https://github.com/ms609/TreeTools/issues/149
tree$node.label <- paste0("Node_", 1:tree$Nnode)

# Note that tips are added only to the dated tree. Also the internal
# node labels are regenerated only for the dated tree. The
# rests of the list elements within the dated tree object remain unchanged.

# export tree
saveRDS(tree, file = "dated_tree.rds")
