rm (list = ls())

args <- commandArgs(trailingOnly = TRUE)

# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

tree_paths <- args[1:length(args)]

trees <- list()
for (i in seq_along(tree_paths)) {
    trees[[i]] <- readRDS(tree_paths[i])
}
names(trees) <- sapply(trees, function(x) x$root_method)

ll <- sapply(trees, function(x) x$loglik)

index <- which(ll == max(ll))[1]

tree <- trees[[index]]

if (grepl("_[0-9]$", tree$root_method)) {
    tree$root_method <- gsub("_[0-9]$", "", tree$root_method)
}

saveRDS(tree, "final_dated_tree.rds")

# end logging
if (!interactive()) {
  sink(con)
}
