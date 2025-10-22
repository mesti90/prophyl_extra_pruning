library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-ts", "--trees_shrink"),
    type = "character",
    help = "A list of dated trees in rds format (from treeshrink)."
  )
  make_option(
    c("-tp", "--trees_prune"),
    type = "character",
    help = "A list of dated trees in rds format (from treepruner)."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    trees_shrink = "dated_trees.rds"
    trees_prune = "dated_trees.rds"
  )
}

# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

trees_shrink <- readRDS(args$trees_shrink)
trees_prune <- readRDS(args$trees_prune)
trees <- c(trees_shrink, trees_prune)

# choose tree with highest log likelihood
ll <- sapply(trees, function(x) x$loglik)
index <- which(ll == max(ll))[1]
tree <- trees[[index]]

if(index <= length(trees_shrink)) {
    method_used <- "treeshrink"
} else {
    method_used <- "treepruner"
}

# choose first tree where rooting approach was RTT RMS
#root_method <- sapply(trees, function(x) x$root_method)
#index <- which(root_method == "rtt_rms_1")
#tree <- trees[[index]]

if (grepl("_[0-9]$", tree$root_method)) {
    tree$root_method <- gsub("_[0-9]$", "", tree$root_method)
}

saveRDS(tree, "final_dated_tree.rds")

writeLines(method_used, "selected_pruning_method.txt")
# end logging
if (!interactive()) {
  sink(con)
}
