#' Add state change counts
#' 
#' This function counts the number of ancestral state changes between a
#' common ancestor and each tip which descends from this common ancestor. The 
#' function currently counts the number of state changes from root.
#' @examples 
#' \dontrun{
#' add_scc(tree_tbl, "k_serotype")
#' }
#' @export
add_scc <- function(tree_tbl, target) {
  index_root <- which(tree_tbl$parent == tree_tbl$node)
  if (length(index_root) == 0) {
    stop("Root node seems to be missing. Check.")
  }
  if (length(index_root) > 1) {
    stop("Multiple root nodes detected. Check.")
  }
  tree <- ape::as.phylo(treeio::as.treedata(tree_tbl))
  scc <- vector()
  for (i in 1:nrow(tree_tbl)) {
    if (i == index_root) scc <- c(scc, 0) else {
      index <- ape::nodepath(tree, from = index_root, to = tree_tbl$node[i])
      series <- tree_tbl[[target]][index]
      scv <- vector()
      for (j in 1:(length(series)-1)) {
        scv <- c(scv, scp(from = series[j], to = series[j+1])$prob)
      }
      scc <- c(scc, sum(scv))
    }
  }
  index_target <- which(names(tree_tbl) == target)
  new_tbl <- dplyr::bind_cols(
    tree_tbl[,1:(index_target+2)],
    scc = scc,
    tree_tbl[,(index_target+3):ncol(tree_tbl)]
  )
  names(new_tbl)[index_target+3] <- paste0(target, "_sc")
  return(new_tbl)
}
