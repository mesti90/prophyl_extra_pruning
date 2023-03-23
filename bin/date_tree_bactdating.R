rm(list = ls())

if (!interactive()) {
  # create log file and start logging
  con <- file("log.txt")
  sink(con, split = TRUE)
}

Sys.time()

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5) {
    stop("Script requires 5 arguments. Check.")
  }
  project_dir <- args[1]
  tree_path <- args[2]
  snps_path <- args[3]
  assemblies_path <- args[4]
  branch_dimension <- args[5]
} else {
  test_dir <- "~/Methods/prophyl-tests/test-date_tree_bactdating"
  project_dir <- "~/Methods/prophyl"
  tree_path <- paste0(test_dir, "/treeshrink.tre")
  snps_path <- paste0(test_dir, "/shrinked_snps.fasta")
  assemblies_path <- paste0(test_dir, "/assemblies.tsv")
  branch_dimension <- "snp_per_genome"
  sample_size <- 100
}

library(ape)
library(BactDating)
library(data.table)
library(devtools)
library(lubridate)
library(tidytree)
library(tidyverse)

load_all(project_dir)

# Input validation
branch_dimension <- match.arg(
  branch_dimension,
  choices = c("snp_per_genome", "snp_per_site")
)

tree <- ape::read.tree(tree_path) %>%
  ape::multi2di() %>%
  ape::unroot()

if(interactive()){
  if (length(tree$tip.label) > sample_size) {
    set.seed(0)
    index <- sample(1:length(tree$tip.label), sample_size)
    tree <- ape::keep.tip(tree, tree$tip.label[index])
  }
}

if (branch_dimension == "snp_per_site") {
  # rescale branch lengths from per site to per genome
  # ?BactDating::bactdate()
  tree$edge.length <- tree$edge.length * alignment_length
}

assemblies <- read.csv(assemblies_path, sep = "\t")

# drop tips which cannot be found in assembly table and give a warning
index <- which(tree$tip.label %in% assemblies$assembly == FALSE)
if (length(index) > 0) {
  tips_to_drop <- tree$tip.label[index]
  tips_to_drop_collapsed <- paste(tips_to_drop, collapse = ", ")
  tree <- ape::drop.tip(tree, tips_to_drop)
  msg <- paste0(
    "One or more tips could not be found in assembly table and were dropped: ",
    tips_to_drop_collapsed,
    "."
  )
  warning(msg)
}

# filter to assemblies that are included in the tree
index <- which(assemblies$assembly %in% tree$tip.label == FALSE)
if (length(index) > 0) {
  assemblies <- assemblies[-index, ]
}

# filter to relevant columns
assemblies <- assemblies[, which(names(assemblies) %in% c(
  "assembly", "collection_date", "collection_day"
))]

# add new temporal variable
assemblies$date <- date_middle(assemblies$collection_date)
assemblies$date <- lubridate::decimal_date(assemblies$date)

# rearrange rows to match tip labels
index <- sapply(tree$tip.label, function(x) {
  which(assemblies$assembly == x)
})
assemblies <- assemblies[index, ]

# TODO add tests for functions add similar checks to scripts
# check rearrangement was successful
if (all(tree$tip.label == assemblies$assembly) == FALSE) {
  stop("Tree tip labels and metadata rows are not in the same order.")
}

dated_tree <- bactdate(
  tree = tree, 
  date = assemblies$date,
  model = "arc"
)

if (!interactive()) {
  # export tree object
  saveRDS(dated_tree, file = "dated_tree.rds")
  # export tree
  ape::write.tree(
    dated_tree,
    file = "dated_tree.tre"
  )
  # export trace
  pdf("trace.pdf", width = 10, height = 6)
  plot(dated_tree, 'trace')
  dev.off()
  # export root-to-tip regression plot
  pdf("root_to_tip_regression.pdf", width = 10, height = 6)
  roottotip(tree, assemblies$date)
  dev.off()
}

if (!interactive()) {
  # end logging
  sink(con)
}
