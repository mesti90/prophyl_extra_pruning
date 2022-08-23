#' Phylogenetic distances
#' 
#' This function calculates phylogenetic (patristic) distances for all pairs of
#' tips within a phylogenetic tree.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @return a tibble with three columns.
#' @examples
#' \dontrun {
#' # TODO put a tree_tbl into data so it can be imported.
#' phylodist(tree_tbl)
#' }
phylodist <- function(tree_tbl) {
  tree <- ape::as.phylo(treeio::as.treedata(tree_tbl))
  phylodist <- ape::cophenetic.phylo(tree)
  phylodist[lower.tri(phylodist)] <- NA
  diag(phylodist) <- NA
  phylodist_df <- dplyr::bind_cols(
    "tip1" = row.names(phylodist),
    as.data.frame(phylodist)
  )
  phylodist_df <- tidyr::pivot_longer(
    phylodist_df,
    cols = 2:ncol(phylodist_df),
    names_to = "tip2",
    values_to = "phylodist"
  )
  phylodist_df <- phylodist_df[which(!is.na(phylodist_df$phylodist)),]
  return(phylodist_df)
}
