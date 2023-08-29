# This script is used to root a subset tree using the same method as the tree
# without subsetting. The script implements a generalised function which
# dispathes rooting to the appropriate function based on the root_method. If
# the root method is rtt, then the function will return the top 5 trees.

rm(list = ls())

# create log file and start logging
con <- file("log.txt")
sink(con, split = TRUE)

library(devtools)
library(dplyr)
library(ggplot2)
library(optparse)

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to the project directory."
  ),
  make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "Path to the assemblies file."
  ),
  make_option(
    c("-d", "--dated_tree"),
    type = "character",
    help = "Path to a dated pylogenetic tree which includes all tips, from treedater, in .rds format, includes the root_method as well."
  ),
  make_option(
    c("-s", "--subset_tree"),
    type = "character",
    help = "Path to a non-dated phlogenetic tree for the subset, in .nwk format."
  ),
  make_option(
    c("-t", "--threads"),
    type = "integer",
    default = 1,
    help = "Number of threads to use."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/aci/prophyl",
    assemblies = "assemblies.tsv",
    dated_tree = "final_dated_tree.rds",
    subset_tree = "subsample_004.nwk",
    threads = 10
  )
}

# load custom functions from the project directory
load_all(args$project_dir)

# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

# read big tree
big_tree <- readRDS(args$dated_tree)
# extract root method, the same method wil be used for all subset trees
root_method <- big_tree$root_method

# read subset tree
tree <- ape::read.tree(args$subset_tree)

# if tree is rooted, unroot
if (ape::is.rooted(tree)) {
  tree <- ape::unroot(tree)
}

# read assemblies
assemblies <- read.csv(args$assemblies, sep = "\t")

# collect tip dates in the same order as tree$tip.label
tip_dates <- unname(sapply(tree$tip.label, function(x) {
  index <- which(assemblies$assembly == x)
  assemblies$collection_day[index]
}))

# convert tip dates to numeric for root_tree()
tip_dates <- as.numeric(as.Date(tip_dates))

# TODO: look for better objectives
objective_rlm_slope <- function(x,y) MASS::rlm(y ~ x)$coef[2]
objective_rlm_rms <- function(x,y) -summary(MASS::rlm(y ~ x))$sigma^2

objective <- list(
  "correlation" = NULL,
  "rsquared" = NULL,
  "rms" = NULL,
  "rlm_slope" = objective_rlm_slope,
  "rlm_rms" = objective_rlm_rms
)

# TODO: make this a generalised function and move it to R directory
root_tree <- function(tree, root_method, tip_dates) {
  if (root_method == "midpoint") {
    rooted_tree <- phytools::midpoint.root(tree)
  } else if (root_method == "mad") {
    rooted_tree <- root_mad(
      tree,
      output_mode = "full",
      cache = TRUE,
      threads = args$threads,
      verbose = TRUE
    )
  } else if (grepl("^rtt", root_method)) {
    root_method <- gsub("^rtt_", "", root_method)
    index <- which(names(objective) == root_method)
    if (length(index) == 0) {
      stop("Objective not found.")
    }
    rooted_tree <- try(root_rtt(
      t = tree,
      tip.dates = tip_dates,
      topx = 5,
      ncpu = args$threads,
      objective = names(objective)[[index]],
      objective_fn = objective[[index]]
    ), silent = TRUE)
    if (inherits(rooted_tree, "try-error")) {
      warning("RTT failed, falling back to midpoint root.")
      rooted_tree <- phytools::midpoint.root(tree)
    }
  } else {
    stop("Root method not found.")
  }
  return(rooted_tree)
}

rooted_tree <- root_tree(tree, root_method, tip_dates = tip_dates)
names(rooted_tree) <- paste0(root_method, "_", 1:5)

subset_id <- strsplit(basename(args$subset_tree), "\\.")[[1]][1]

# export rooted tree object
saveRDS(rooted_tree, file = paste0("rooted_trees_", subset_id, ".rds"))

# end logging
if (!interactive()) {
  sink(con)
}
