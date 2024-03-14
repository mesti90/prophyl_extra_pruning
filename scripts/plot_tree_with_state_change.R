library(ggnewscale)
library(ggimage)
library(ggplot2)
library(ggtree)
library(optparse)
library(qualpalr)
rm(list = ls())

# dates from Hun to Eng format
Sys.setlocale("LC_TIME", "C")

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "The project directory."
  ),
  make_option(
    c("-t", "--tree_tbl"),
    type = "character",
    help = "A tree in table format with predicted ancestral states."
  ),
  make_option(
    c("-T", "--target"),
    type = "character",
    help = "A column name in the tree table to be used as the target variable."
  ),
  make_option(
    c("-h", "--heatmap_vars"),
    type = "character",
    help = "A comma separated list of variables that should be used for
    generating heatmaps. Currently only a single variable is supported."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci/prophyl",
    tree_tbl = "data/global_ST2_tree/tree_tbl.rds",
    target = "city_pooled",
    heatmap_vars = "city_pooled"
  )
}

library(devtools)
load_all(args$project_dir)

tree_tbl_file <- args$tree_tbl

tree_tbl <- readRDS(tree_tbl_file)

index <- which(is.na(tree_tbl$branch.length))

if (length(index) != 1) {
  stop("Could not find tree root.")
}

set.seed(0)
geo_cols <- get_colors(tree_tbl[1:(index-1), ], args$target)

tree_tbl <- dplyr::left_join(tree_tbl, geo_cols, by = args$target)

testthat::expect_equal(sum(is.na(tree_tbl$color[1:(index-1)])), 0)

tree <-  treeio::as.treedata(tree_tbl)

max_date <- max(tree_tbl$collection_day, na.rm = TRUE)

mat <- as.matrix(tree_tbl[, args$heatmap_vars])
rownames(mat) <- tree_tbl$label
colnames(mat) <- args$heatmap_vars

index_ambiguous_nodes <- grep("\\|", tree_tbl[[args$target]])

index_transmission_nodes <- vector()
for (i in 1:nrow(tree_tbl)) {
  target_child <- tree_tbl[[args$target]][i]
  target_parent <- tree_tbl[[args$target]][which(tree_tbl$node == tree_tbl$parent[i])]
  if (target_child == target_parent) next() else {
    index_transmission_nodes <- c(index_transmission_nodes, i)
  }
}

p <- ggtree(tree, aes(color = get(args$target)), mrsd = max_date) +
  theme_tree2()+
  scale_x_ggtree()+
  scale_color_manual(
    values = geo_cols$color,
    limits = geo_cols[[args$target]]
  )+
  geom_point2(
    aes(subset = (node %in% index_ambiguous_nodes)),
    shape = 21,
    size = 10,
    fill = "orange"
  )+
  geom_label2(
    aes(
      x = branch,
      subset = (node %in% index_transmission_nodes),
      label = "INTRO?"
    ),
    size = 3,
    fill = "red")+
  geom_tiplab(align = TRUE, size = 2)+
  geom_label(aes(label = city_pooled), size = 2)+
  theme(legend.position = "none")

cls <- geo_cols$color
names(cls) <- geo_cols$city_pooled

p2 <- gheatmap(
  p,
  mat,
  offset = 1,
  width = 0.02,
  colnames_offset_y = 0,
  colnames_position = "top",
  legend_title = args$heatmap_vars
  )+
  labs(fill = args$heatmap_vars)+
  scale_fill_manual(values = cls)+
  guides(colour = "none", fill = "none")

ggsave(
  file = "tree_with_state_changes.pdf",
  height = 0.2*nrow(tree_tbl),
  width = 0.01*nrow(tree_tbl),
  units = "cm",
  limitsize = FALSE
)

#ggsave(
#  file = "tree_with_state_changes.png",
#  height = 8*nrow(tree_tbl),
#  width = 4*nrow(tree_tbl)
#)
