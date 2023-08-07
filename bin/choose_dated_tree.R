library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-t", "--trees"),
    type = "character",
    help = "A list of dated trees in rds format."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    trees = "dated_trees.rds"
  )
}

# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

trees <- readRDS(args$trees)

# choose tree with highest log likelihood
# ll <- sapply(trees, function(x) x$loglik)
# index <- which(ll == max(ll))[1]
# tree <- trees[[index]]

# choose first tree where rooting approach was RTT RMS
root_method <- sapply(trees, function(x) x$root_method)
index <- which(root_method == "rtt_rms_1")
tree <- trees[[index]]

if (grepl("_[0-9]$", tree$root_method)) {
    tree$root_method <- gsub("_[0-9]$", "", tree$root_method)
}

saveRDS(tree, "final_dated_tree.rds")

# end logging
if (!interactive()) {
  sink(con)
}
