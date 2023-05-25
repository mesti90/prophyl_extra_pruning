#' Remove Tips in a Phylogenetic Tree
#' 
#' This function is a wrapper around \code{ape::drop.tip()} to facilitate
#' dropping tips from a phylogenetic tree in tibble format.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param tip a vector of mode numeric or character specifying the tips to delete
#' @param ... additional arguments passed on to \code{ape::drop.tip()}.
#' @examples 
#' \dontrun {
#' small_tbl <- drop_tip(tree_tbl, tips = c("GCA_016487125.1","GCA_016487085.1"))
#' }
#' @export
drop_tip <- function(tree_tbl, tip, ...) {
  tree <- ape::as.phylo(treeio::as.treedata(tree_tbl))
  small <- ape::drop.tip(tree, tip = tip, ...)
  small_tbl <- tibble::as_tibble(small)
  small_tbl <- dplyr::left_join(
    small_tbl,
    tree_tbl[,-which(names(tree_tbl) %in% c("parent", "node", "branch.length"))],
    by = "label"
  )
  return(small_tbl)
}
# NOTE: not sure if the function does the right thing, check.