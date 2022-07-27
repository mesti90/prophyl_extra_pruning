rm(list=ls())
library(ggnewscale)
library(ggimage)
library(ggplot2)
library(ggtree)
library(qualpalr)

# dates from Hun to Eng format
Sys.setlocale("LC_TIME", "C")

args <- commandArgs(trailingOnly = TRUE)
tree_tbl_file <- args[1]

tree_tbl <- readRDS(tree_tbl_file)

# instead of country make this more generic

country <- tree_tbl$country[which(!is.na(tree_tbl$country))]

country_cols <- data.frame(
  country = sort(unique(country)),
  country_col = qualpalr::qualpal(length(unique(country)), colorspace = "pretty_dark")$hex
)

tree_tbl <- dplyr::left_join(tree_tbl, country_cols, by = "country")
tree_tbl$country_col[which(is.na(tree_tbl$country_col))] <- "#D3D3D3"

ccode <- function(country) {
  if (country == "kosovo") return("XK")
  if (country == "south_korea") return("KR")
  countrycode::countrycode(country, "country.name", "iso2c")
}

tree_tbl$country_short <- unname(sapply(tree_tbl$country, function(x) {
  candidates <- strsplit(x, split = "\\|")[[1]]
  shorts <- vector()
  for (i in candidates) {
    shorts <- c(
      shorts, 
      paste(ccode(i), collapse = "")
    )
  }
  paste(shorts, collapse = "|")
}))


tree <-  treeio::as.treedata(tree_tbl)

max_date <- max(tree_tbl$collection_day, na.rm = TRUE)

mat <- as.matrix(tree_tbl[,"k_serotype"])
rownames(mat) <- tree_tbl$label
colnames(mat) <- "k_serotype"

index_ambiguous_nodes <- grep("\\|", tree_tbl$country)

index_transmission_nodes <- vector()
for (i in 1:nrow(tree_tbl)) {
  country_child <- tree_tbl$country[i]
  country_parent <- tree_tbl$country[which(tree_tbl$node == tree_tbl$parent[i])]
  if (country_child == country_parent) next() else {
    index_transmission_nodes <- c(index_transmission_nodes, i)
  }
}

country_cols <- tree_tbl[,which(names(tree_tbl) %in% c("country", "country_col"))]
country_cols <- dplyr::distinct(country_cols)
country_cols$country <- unname(country_cols$country)

p <- ggtree(tree, aes(color = country), mrsd = max_date) +
  theme_tree2()+
  scale_x_ggtree()+
  scale_color_manual(
    values = country_cols$country_col,
    limits = country_cols$country
  )+
  geom_point2(
    aes(subset = (node %in% index_ambiguous_nodes)),
    shape = 21,
    size = 20,
    fill = "orange"
  )+
  geom_label2(
    aes(
      x = branch,
      subset = (node %in% index_transmission_nodes),
      label = "INTRO?"
    ),
    size = 5,
    fill = "red")+
  geom_tiplab(align = TRUE)+
  geom_label(aes(label = country_short))+
  theme(legend.position = "none")

cls <- qualpal(length(unique(tree_tbl$k_serotype)), colorspace = "pretty")$hex
p2 <- gheatmap(
  p,
  mat,
  offset = 1,
  width = 0.02,
  colnames_offset_y = 0,
  colnames_position = "top",
  legend_title = "K serotype"
  )+
  labs(fill = "k_serotype")+
  scale_fill_manual(values = cls)+
  guides(colour = "none")

pdf(
  file = "dated_tree.pdf",
  height = 0.1*nrow(tree_tbl),
  width = 0.05*nrow(tree_tbl),
  compress = FALSE
)
p2
dev.off()

png(
  file = "dated_tree.png",
  height = 8*nrow(tree_tbl),
  width = 4*nrow(tree_tbl)
)
p2
dev.off()
