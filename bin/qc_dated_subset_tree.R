library(dplyr)
library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("-s", "--subsample_id"),
    type = "character",
    help = "Subset ID"
  ),
  make_option(
    c("-t", "--tree"),
    type = "character",
    help = "A dated tree in rds format."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    subset_id = "subsample_001",
    project_dir = "~/Methods/prophyl",
    tree = "dated_tree.rds"
  )
}

library(devtools)
load_all(args$project_dir)

tree <- readRDS(args$tree)

df <- data.frame(
    subsample_id = args$subsample_id,
    mrca = tree$timeOfMRCA
)

write.table(
    df,
    file = paste0(args$subsample_id, "_qc.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
