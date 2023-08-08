# This script is one of multiple scripts that plot phylogenetic trees. This
# script creates a phylogenetic tree plot with fan layout. The tree will also 
# contain heatmaps for MLST, K type, country and city and will highlight tips
# that were included in phage-host laboratory tests.
#
# The script takes
# - a dated phylogenetic tree which was created by the Nextflow pipeline and
#   which includes any duplicates that were removed during tree building
# - a metadata table which contains typing results and geographical locations
#   for aci strains
# - an rds file which contains a heatmap abour phage-host laboratory test
#   results
# 
# It is expected that this script will be run after the Nextflow pipeline which
# creates the dated phylogenetic trees. The working directory for running this
# script should be the directory from which the Nextflow pipeline that build the
# phylogenetic tree was executed.
#
# The dated phylogenetic tree will be automatically created by the Nextflow
# pipeline and will be located in the results/add_duplicates directory. The rest
# of the files must be supplied manually, i.e. they must be present in the
# working directory.

rm(list = ls())

library(devtools)
library(dplyr)
library(treeio)

load_all("~/Methods/prophyl")

# path to the phylogenetic tree
tree_path <- "./results/add_duplicates/dated_tree.rds"
tree <- readRDS(tree_path)

# path to the file used for creating the phylogenetic tree. This file also
# contains metadata required by this script.
meta_path <- "assemblies.tsv"
meta <- read.csv("assemblies.tsv", sep = "\t")

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

# path to the file which contains assemblies that should be highlighted
labres_path <- "heatmap_phages_PFU.rds"
labres <- readRDS(labres_path)

# filter to unique assemblies
labres <- labres$data[,c("Assembly", "MLST", "KL")]
labres <- dplyr::distinct(labres)

ST <- unique(meta$mlst)

if (length(ST) > 1) {
  stop("Script is designed to work with a single ST.")
}

labres <- labres[which(labres$MLST == ST),]

# Throw a message if tested assemblies are missing from the tree

# missing assemblies
missing_assemblies <- labres$Assembly[which(
  labres$Assembly %in% meta$assembly == FALSE
)]

# print informative message
if (length(missing_assemblies) > 0) {
  missing_assemblies_collapsed <- paste(
    missing_assemblies, collapse = ", "
  )
  msg <- paste0(
    "Some assemlies are missing from the tree: ",
    missing_assemblies_collapsed, "."
  )
  message(msg)
}

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

k_types <- unique(labres$KL)

# prepare a tree for each K type
for (k in k_types) {
  # subset to those with given K type
  # subset to those from Europe
  # subset to those that are on the tree
  meta_small <- meta[which(
    meta$k_serotype == k &
    meta$continent == "europe" &
    meta$assembly %in% tree$tip.label
  ),]
  small_tree <- ape::keep.tip(tree, meta_small$assembly)
  small_tree_tbl <- treeio::as_tibble(small_tree)
  small_tree_tbl <- dplyr::left_join(
    small_tree_tbl,
    meta_small[,c(
      "assembly", "continent", "country", "city", "mlst",
      "k_serotype"
    )],
    by = c("label" = "assembly")
  )
  small_tree_tbl$tested <- FALSE
  small_tree_tbl$tested[which(small_tree_tbl$label %in% labres$Assembly)] <- TRUE
  s1 <- 0
  s2 <- 1
  while (s2 > s1) {
    s1 <- sum(small_tree_tbl$tested, na.rm = TRUE)
    parents <- unique(small_tree_tbl$parent[which(small_tree_tbl$tested == TRUE)])
    index <- which(small_tree_tbl$node %in% parents)
    small_tree_tbl$tested[small_tree_tbl$node %in% index] <- TRUE
    s2 <- sum(small_tree_tbl$tested, na.rm = TRUE)
  }
  # use single color for MLST
  mlst_colors <- data.frame(
    mlst = ST,
    color = "grey50"
  )
  # use single color for K type
  k_colors <- data.frame(
    k_serotype = k,
    color = "grey50"
  )
  # define colors for countries
  # TODO use defined country colors where possible
  country_colors <- data.frame(
    country = sort(unique(meta_small$country)),
    color = qualpalr::qualpal(length(unique(meta_small$country)), "pretty")$hex
  )
  # define colors for cities
  # only include cities which come from countries with multiple cities
  citycount <- small_tree_tbl %>% 
    group_by(country) %>% 
    summarise(count = length(unique(city[which(!is.na(city))])))
  selected_countries <- citycount$country[which(citycount$count > 1)]
  if (length(selected_countries) == 0) {
    index <- numeric()
  } else {
    index <- which(small_tree_tbl$country %in% selected_countries)
  }
  if (length(index) == 0) {
  } else {
    selected_cities <- unique(small_tree_tbl$city[index])
    selected_cities <- selected_cities[which(!is.na(selected_cities))]
    selected_cities_colors <- qualpalr::qualpal(length(selected_cities), "pretty")$hex
  }
  city_colors <- data.frame(
    city = sort(selected_cities),
    color = selected_cities_colors
  )
  # for ST492 use
  # - heatmap_width = 0.75
  # for all other STs, use
  # - heatmap_width = 2.5
  g <- plot_tree_fan(
    small_tree_tbl,
    highlight_var = "tested",
    open_angle = 15,
    heatmap_var = c("mlst", "k_serotype", "country", "city"),
    heatmap_colors = list(
      "mlst" = mlst_colors,
      "k_serotype" = k_colors,
      "country" = country_colors,
      "city" = city_colors
    ),
    heatmap_width = 2.5,
    heatmap_colnames_font_size = 3,
    verbose = TRUE
  )
  #g1 <- g + guides(
  #  country = guide_legend(order = 1, ncol = 1),
  #  city = guide_legend(order = 2, ncol = 2),
  #  mlst = guide_legend(order = 3),
  #  k_serotype = guide_legend(order = 4),
  #  tested = guide_legend(order = 5)
  #)
  file_name <- paste0(ST, "_", k, "_Europe.pdf")
  ggsave(
    filename = file_name,
    limitsize = FALSE,
    width = 32,
    height = 20,
    units = "cm"
  )
}
