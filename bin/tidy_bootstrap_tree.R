rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

library(ggtree)

tree <- ape::read.tree(args[1])

tree_tbl <- tibble::as_tibble(tree)
tree_tbl <- tree_tbl[, -which(names(tree_tbl) == "branch.length")]
tree_tbl$bs <- suppressWarnings(round(as.numeric(tree_tbl$label)/100, 2))
tree_tbl$label <- ifelse(is.na(tree_tbl$bs), tree_tbl$label, NA)

tree_tbl$label[which(tree_tbl$parent == tree_tbl$node)] <- "Node_1"

class(tree_tbl) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")

saveRDS(tree_tbl, file = "tree_tbl.rds")
