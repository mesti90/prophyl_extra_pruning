#' Plot a phlyogenetic tree
#' 
#' This function is used to plot a phylogenetic tree together with predicted
#' ancestral state changes that may have occurred along the branches of the
#' tree.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param state_var character; the name of the variable which stores the
#' character state.
#' @param date_var character; the name of the variable which stores sample
#' collection dates.
#' @param node_lab charcater; the name of the variable which stores node labels.
#' If \code{NULL}, node labels will be character states according to
#' \code{state_var}.
#' @param export logical; should the plot be exported?
#' @import ggnewscale
#' @import ggimage
#' @import ggplot2
#' @import ggtree
#' @importFrom qualpalr qualpal
#' @export
plot_tree <- function(tree_tbl,
                      state_var,
                      date_var = "collection_day",
                      node_lab = NULL,
                      export = TRUE) {
  tree <-  treeio::as.treedata(tree_tbl)
  max_date <- max(tree_tbl$collection_day, na.rm = TRUE)
  cls <- qualpalr::qualpal(
    length(unique(tree_tbl[[state_var]])), colorspace = "pretty")$hex
  index_ambiguous_nodes <- grep("\\|", tree_tbl[[state_var]])
  index_transmission_nodes <- vector()
  for (i in 1:nrow(tree_tbl)) {
    child <- tree_tbl[[state_var]][i]
    parent <- tree_tbl[[state_var]][which(tree_tbl$node == tree_tbl$parent[i])]
    if (child == parent) next() else {
      index_transmission_nodes <- c(index_transmission_nodes, i)
    }
  }
  if (is.null(node_lab)) {
    lab <-  state_var
  } else {
    lab <- node_lab
  }
  mat <- as.matrix(tree_tbl[[state_var]])
  rownames(mat) <- tree_tbl$label
  colnames(mat) <- state_var
  p <- ggtree(tree, aes(color = get(state_var)), mrsd = max_date) +
    #theme_tree2()+
    #scale_x_ggtree()+
    scale_color_manual(values = cls, limits = unique(tree_tbl[[state_var]]))+
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
        label = "CHANGE"
      ),
      size = 5,
      fill = "red")+
    geom_tiplab(align = TRUE)+
    geom_label(aes(label = get(lab)))+
    theme(legend.position = "none")
  p2 <- gheatmap(
    p,
    mat,
    offset = 1.2,
    width = 0.02,
    colnames_offset_y = 0,
    colnames_position = "top",
    legend_title = state_var
    )+
    labs(fill = state_var)+
    scale_fill_manual(values = cls, limits = unique(tree_tbl[[state_var]]))+
    guides(colour = "none")
  
  if (export == TRUE) {
    ggsave(
      filename = "tree.pdf",
      height = 20*nrow(tree_tbl),
      width = 8*nrow(tree_tbl),
      units = "px",
      limitsize = FALSE
    )
    message("Plot exported as tree.pdf.")
  } else {
    p2
  }
}
