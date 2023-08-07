library(dplyr)
library(optparse)
library(TreeTools)
rm(list = ls())

args_list <- list(
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
    tree = "final_dated_tree.rds",
    duplicates = "duplicates.txt"
  )
}

# import tree
tree <- readRDS(args$tree)

# import duplicates
duplicates <- readLines(args$duplicates)

if (length(duplicates) > 0) {
  
  # format duplicates
  duplicates <- duplicates %>%
    gsub("^[0-9]+\\\t", "", .) %>%
    strsplit(., ", ")
  names(duplicates) <- sapply(duplicates, function(x) x[1])
  
  # add duplicates to tree
  for (i in duplicates) {
    for (j in 2:length(i)) {
      tree <- TreeTools::AddTip(
        tree,
        where = i[1],
        label = i[j],
        edgeLength = 0
      )
    }
  }
  
  # consistency check
  duplicate_tips <- duplicates %>% unlist() %>% unique()
  testthat::expect_true(all(duplicate_tips %in% tree$tip.label))
  
}

# export tree
saveRDS(tree, file = "dated_tree.rds")
