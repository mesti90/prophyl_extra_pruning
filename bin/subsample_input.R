# This script takes a table of assemblies and draws samples using the selected
# sampling strategy. Each subsampled table will then be used to build and date a
# subset tree. Each subset tree will be used to calculate a number of data
# points for the relative risk plots.

library(optparse)
rm(list=ls())

args_list <- list(
  make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "Path to assemblies file."
  ),
  make_option(
    c("-t", "--tree"),
    type = "character",
    help = "A dated tree in rds format."
  ),
  make_option(
    c("-c", "--subsample_count"),
    type = "character",
    help = "Number of subsample sets to draw."
  ),
  make_option(
    c("-C", "--subsample_tipcount"),
    type = "character",
    help = "Number of tips to draw in each subsample set."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    assemblies = "assemblies.tsv",
    tree = "dated_tree.rds",
    subsample_count = 10,
    subsample_tipcount = 10
  )
}

# read assemblies
assemblies <- read.csv(args$assemblies, sep = "\t")
# read tree
tree <- readRDS(args$tree)
# number of subsample sets to draw
subsample_count <- as.numeric(args$subsample_count)
# number of tips to draw in each subsample set
subsample_tipcount <- as.numeric(args$subsample_tipcount)

# The shrinked tree may contain less tips than the original tree
# Only sample assemblies that are included in the shrinked tree
index <- which(assemblies$assembly %in% tree$tip.label == FALSE)
if (length(index) > 0) {
  assemblies <- assemblies[-index, ]
}

# sampling strategy
# can be either "random", "balanced" or "focused".
type <- "random"
# a variable within the input table used for balancing.
# only used if the sampling strategy is "balanced"
balance_by <- NULL
# a variable within the input table used for focusing.
# only used if sampling strategy is "focused".
focus_by <- "continent"
focus_on <- "europe"
# ratio of samples to take from focus group.
focus_ratio <- 0.75

if (!is.null(balance_by) && balance_by %in% names(assemblies) == FALSE) {
  stop("Subsampling failed, variable '", balance_by, "' not found.")
}

if (type == "random") {
  set.seed(0)
  digits <- ceiling(log10(subsample_count+1))
  for (i in 1:subsample_count) {
    zeroes <- digits - floor(log10(i))
    filename = paste0(c("subsample_", rep(0, times = zeroes), i, ".tsv"), collapse = "")
    write.table(
      assemblies[sample(1:nrow(assemblies), subsample_tipcount, replace = FALSE), ],
      file = filename,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}

if (type == "balanced") {
  set.seed(0)
  years_min <- min(assemblies$collection_year, na.rm = TRUE)
  years_max <- max(assemblies$collection_year, na.rm = TRUE)
  years <- years_min:years_max
  
  size = 0
  k = 0
  
  while (size < subsample_tipcount & size < nrow(assemblies)) {
    k = k + 1
    index <- vector()
    for (i in years) {
      bins <- unique(assemblies[[balance_by]][which(assemblies$collection_year == i)])
      for (j in bins) {
        idx <- which(assemblies$collection_year == i & assemblies[[balance_by]] == j)
        if (length(idx) >= k) {
          index <- c(index, sample(idx, k))
        } else {
          index <- c(index, idx)
        }
      }
    }
    size <- length(index)
  }
  if (k == nrow(assemblies)) {
    stop("Subsampling did not converge. Check.")
  } else {
    digits <- ceiling(log10(subsample_count+1))
    for (h in 1:subsample_count) {
      index <- vector()
      for (i in years) {
        bins <- unique(assemblies[[balance_by]][which(assemblies$collection_year == i)])
        for (j in bins) {
          idx <- which(assemblies$collection_year == i & assemblies[[balance_by]] == j)
          if (length(idx) >= k) {
            index <- c(index, sample(idx, k))
          } else {
            index <- c(index, idx)
          }
        }
      }
      zeroes <- digits - floor(log10(h))
      filename = paste0(c("subsample_", rep(0, times = zeroes), h, ".tsv"), collapse = "")
      write.table(
        assemblies[index, ],
        file = filename,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  }
}

if (type == "focused") {
  set.seed(0)
  focus_count <- round(focus_ratio * subsample_tipcount, 0)
  no_focus_count <- subsample_tipcount - focus_count
  if (focus_count == 0 | no_focus_count == 0) {
    stop("Number of assemblies in and outside of focus group must be non zero.")
  }
  digits <- ceiling(log10(subsample_count + 1))
  for (i in 1:subsample_count) {
    zeroes <- digits - floor(log10(i))
    filename = paste0(
      c("subsample_", rep(0, times = zeroes), i, ".tsv"), collapse = "")
    focus_index <- which(assemblies[[focus_by]] == focus_on)
    if (length(focus_index) >= focus_count) {
      focus <- sample(
      focus_index,
      focus_count,
      replace = FALSE
      )
    } else {
      stop("Not enough assemblies in focus group to subsample")
    }
    # exclude any assemblies where focus variable is NA
    no_focus_index <- which(assemblies[[focus_by]] != focus_on)
    if (length(no_focus_index) >= no_focus_count) {
      no_focus <- sample(
      no_focus_index,
      no_focus_count,
      replace = FALSE
      )
    } else {
      stop("Not enough assemblies in non-focus group to subsample")
    }
    write.table(
      assemblies[c(focus, no_focus), ],
      file = filename,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}
