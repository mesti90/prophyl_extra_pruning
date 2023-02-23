# What this script does:
# - Plots a phylogenetic tree using the fan layout
# - Tree is filtered to ST2, KL3, Europe
# - The tree may be the whole tree or a random subsample
# - Adds three heatmaps: mlst, k_serotype, country 

rm(list = ls())

library(devtools)
library(dplyr)
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

# filter to ST2, KL3, Europe
meta <- meta[which(meta$mlst == "ST2" & 
                   meta$k_serotype == "KL3" &
                   meta$continent == "europe"),]

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

# Capitalise the first letters of country names
meta$country <- stringr::str_to_title(meta$country)
# Edit some country names manually
meta$country <- gsub("Bosnia_and_herzegovina", "Bosnia and Herzegovina", meta$country)

# Capitalise the first letters of city names
meta$city <- stringr::str_to_title(meta$city)
# Edit some city names manually
meta$city <- gsub("Győr", "Gyor", meta$city)
meta$city <- gsub("Targu_mures", "Targu Mures", meta$city)
meta$city <- gsub("Banja_luka", "Banja Luka", meta$city)

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
if (length(unique(meta$k_serotype)) >= top_k_count) {
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

# define colors for countries
country_colors <- data.frame(
  country = sort(unique(meta$country)),
  color = qualpalr::qualpal(length(unique(meta$country)), "pretty")$hex
)

# define colors for cities
# only include cities which come from countries with multiple cities
citycount <- tree_tbl %>% 
               group_by(country) %>% 
               summarise(count = length(unique(city[which(!is.na(city))])))

selected_countries <- citycount$country[which(citycount$count > 1)]

if (length(selected_countries) == 0) {
  index <- numeric()
} else {
  index <- which(tree_tbl$country %in% selected_countries)
}

if (length(index) == 0) {

} else {
  selected_cities <- unique(tree_tbl$city[index])
  selected_cities <- selected_cities[which(!is.na(selected_cities))]
  selected_cities_colors <- qualpalr::qualpal(length(selected_cities), "pretty")$hex
}

city_colors <- data.frame(
  city = sort(selected_cities),
  color = selected_cities_colors
)

# Exclude mlst and k_serotype from legend
g <- plot_tree_fan(
  tree_tbl,
  drop_tip = tips_to_drop,
  open_angle = 15,
  heatmap_var = c("mlst", "k_serotype", "country", "city"),
  heatmap_colors = list("mlst" = mlst_colors,
                        "k_serotype" = k_colors,
                        "country" = country_colors,
                        "city" = city_colors),
  heatmap_width = 2.5,
  heatmap_colnames_font_size = 3,
  legend_show = c("country", "city"),
  file_name = NULL,
  verbose = TRUE
)

g1 <- g + guides(
  country = guide_legend(ncol = 1),
  city = guide_legend(ncol = 2)
)

ggsave(
      filename = file_name,
      limitsize = FALSE,
      width = 32,
      height = 20,
      units = "cm"
)

#Include mlst and k_serotype in legend
g <- plot_tree_fan(
  tree_tbl,
  drop_tip = tips_to_drop,
  open_angle = 15,
  heatmap_var = c("mlst", "k_serotype", "country", "city"),
  heatmap_colors = list("mlst" = mlst_colors,
                        "k_serotype" = k_colors,
                        "country" = country_colors,
                        "city" = city_colors),
  heatmap_width = 2.5,
  heatmap_colnames_font_size = 3,
  file_name = NULL,
  verbose = TRUE
)

g1 <- g + guides(
  country = guide_legend(ncol = 1),
  city = guide_legend(ncol = 2)
)

file_name2 <- gsub("\\.pdf", "_2.pdf", file_name)

ggsave(
      filename = file_name2,
      limitsize = FALSE,
      width = 32,
      height = 20,
      units = "cm"
)
