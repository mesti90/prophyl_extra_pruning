rm(list = ls())

library(devtools)
library(ggnewscale)
library(ggimage)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(treeio)

args <- commandArgs(trailingOnly = TRUE)

# inputs
# project directory, required for loading functions
project_dir <- args[1]
# a phylogenetic tree in newick format
intree_path <- args[2]
# a metadata file with .rds extension e.g. aci_all.rds
meta_path <- args[3]
# poppunk clusters with .csv extension
pp_path <- args[4]
# sample size, for testing purposes. Will only do sampling if sample size > 0.
sample_size <- args[5]
# output file name
file_name <- args[6]

# parameters
# number of most frequent mlst-s to plot separately, pool the rest as "Other"
top_mlst_count = 6
# number of most frequent k serotypes to plot separately, pool the rest as "Other"
top_k_count = 10
# drop these tips (e.g. becasue they would distort the final plot)
tips_to_drop = c("GCF_014171935.1", "GCA_900495195.1")

# notes
# 1. legend definition based on frequencies 2. (potential) subsampling. Because
# of this the legend will represent most frequent types etc. based on all the 
# data not the subsample.

# script
load_all(project_dir)

tree <- ape::read.tree(intree_path)
meta <- readRDS(meta_path)

# drop tips which do not have related metadata
if (all(tree$tip.label %in% meta$assembly) == FALSE) {
  labels <- tree$tip.label[which(tree$tip.label %in% meta$assembly == FALSE)]
  msg <- paste(labels, collapse = ", ")
  warning("One or more tip labels cannot be found in metadata table: ", msg)
  # drop tips that are not in the table
  tree <- ape::drop.tip(tree, tip = labels)
}

# only keep metadata for assemblies that are on the tree
meta <- meta[which(meta$assembly %in% tree$tip.label),]

# add poppunk clusters
pp <- read.csv(pp_path)
pp <- dplyr::rename(pp, assembly = Taxon, pp = Cluster)

pp$assembly <- gsub("_", ".", pp$assembly)
pp$assembly <- gsub("GCA\\.", "GCA_", pp$assembly)
pp$assembly <- gsub("GCF\\.", "GCF_", pp$assembly)

meta <- dplyr::left_join(meta, pp, by = "assembly")

# Capitalise the first letters of continent names
meta$continent <- stringr::str_to_title(meta$continent)

# collapse mlst and define colors
if (length(unique(meta$mlst)) >= top_mlst_count) {
  top_mlst <- names(sort(table(meta$mlst), decreasing = TRUE))[1:top_mlst_count]
  meta$mlst <- ifelse(meta$mlst %in% top_mlst, meta$mlst, "Other")
  
  mlst_colors <- data.frame(
    mlst = c(sort(top_mlst), "Other"),
    color = c(qualpalr::qualpal(top_mlst_count, "pretty")$hex, "grey50")
  )
} else {
  mlst_colors <- data.frame(
    mlst = sort(unique(meta$mlst)),
    color = ifelse(
      length(unique(meta$mlst)) == 1,
      "grey50",
      qualpalr::qualpal(length(unique(meta$mlst)), "pretty")$hex
    )
  )
}

# collapse k_serotype and define colors
if (length(unique(meta$k_serotype) >= top_k_count)) {
  top_k <- names(sort(table(meta$k_serotype), decreasing = TRUE))[1:top_k_count]
  meta$k_serotype <- ifelse(meta$k_serotype %in% top_k, meta$k_serotype, "Other")
  
  k_colors <- data.frame(
    k_serotype = c(sort(top_k), "Other"),
    color = c(qualpalr::qualpal(top_k_count, "pretty")$hex, "grey50")
  )
} else {
  k_colors <- data.frame(
    k_serotype = sort(unique(meta$k_serotype)),
    color = ifelse(
      length(unique(meta$k_serotype)) == 1,
      "grey50",
      qualpalr::qualpal(length(unique(meta$k_serotype)), "pretty")$hex
    )
  )
}

# define colors for continents

continent_colors <- data.frame(
  continent = sort(unique(meta$continent)),
  color = qualpalr::qualpal(length(unique(meta$continent)), "pretty")$hex
)

# subsample tree
if (sample_size > 0) {
  set.seed(0)
  sample_size <- round(as.numeric(sample_size),0)
  tree <- ape::keep.tip(tree, sample(tree$tip.label, sample_size))
}

tree_tbl <- as_tibble(tree)

tree_tbl <- dplyr::left_join(
  tree_tbl,
  meta[,c(
    "assembly", "continent", "region23", "country", "city", "mlst",
    "k_serotype", "k_confidence", "pp", "xdr"
  )],
  by = c("label" = "assembly")
)

save(list = c("tree_tbl", "tips_to_drop", "mlst_colors", "k_colors", "continent_colors"), file = "test.rda")

plot_tree_fan(
  tree_tbl,
  drop_tip = tips_to_drop,
  open_angle = 15,
  heatmap_var = c("mlst", "k_serotype", "continent"),
  heatmap_colors = list("mlst" = mlst_colors,
                        "k_serotype" = k_colors,
                        "continent" = continent_colors),
  heatmap_colnames_font_size = 3,
  file_name = file_name,
  verbose = TRUE
)
