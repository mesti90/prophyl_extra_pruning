#' Prepare ancestral state network
#' 
#' This function prepares an adjacency matrix which represents a directed
#' ancestral state network along a phylogenetic tree. Nodes in this network
#' represent unique ancestral states, edges represent state changes between
#' ancestors and their descendants.
#' @param tree_tbl tibble; dated tree with predicted ancestral states. This
#' tibble is automatically returned by the nextlow pipeline.
#' @param target character; the target ancestral state to plot. Must match a
#' column name in \code{tree_tbl}.
#' @param cache_file character; the name of the cache file without the file
#' extension. If \code{NULL}, results are not cached.
#' @return an adjacency matrix
#' @seealso plot_snet
#' @examples 
#' \dontrun{
#' snet(tree_tbl, "k_serotype", cache_file = "adjmat_k_serotype")
#' }
#' @export
snet <- function(tree_tbl, target, cache_file = NULL) {
  target <- match.arg(target, names(tree_tbl))
  states <- unique(unname(unlist(strsplit(tree_tbl[[target]], split = "\\|"))))
  if(!is.null(cache_file)) {
    cfpath <- paste0("cache/", cache_file, ".rds")
    if (file.exists(cfpath)) {
      adjmat <- readRDS(cfpath)
      if (mean(rownames(adjmat) == states) != 1) {
        stop("Ancestral states in query and cache do not match. Check.")
      }
      if (mean(colnames(adjmat) == states) != 1) {
        stop("Ancestral states in query and cache do not match. Check.")
      }
    } else {
      adjmat <- matrix(NA, nrow = length(states), ncol = length(states))
      rownames(adjmat) <- states
      colnames(adjmat) <- states
      for (i in 1:length(states)) {
        for (j in 1:length(states)) {
          adjmat[i,j] <- scc(
            states[i], states[j], target = target, tree_tbl = tree_tbl)$count
        }
      }
      if (!dir.exists("cache")) dir.create("cache")
      cfpath <- paste0("cache/", cache_file, ".rds")
      saveRDS(adjmat, file = cfpath)
    }
  }
  return(adjmat)
}
