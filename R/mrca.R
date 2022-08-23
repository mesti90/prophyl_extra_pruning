#' Calculate times to most recent common ancestors
#' 
#' This function calculates phylogenetic (patristic) distances for all pairs of
#' tips within a phylogenetic tree.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param method character; the method to be used for calculating time to the
#' most recent common ancestors. See details for more information.
#' @param ... additional arguments used with \code{colldist()}.
#' @return a tibble with five columns.
#' @details The function implements three methods for calculating times to the
#' most recent common ancestors. \code{"first"} calculates times between the
#' ancestor and the earlier tip, \code("last") calculates times between the
#' ancestor and the more recent tip, \code("average") calculates the average of
#' the two which is identical to the half of the patristic distance.
#' @examples
#' \dontrun {
#' # TODO put a tree_tbl into data so it can be imported.
#' mrca(tree_tbl)
#' }
#' @export
mrca <- function(tree_tbl, method = "average", ...) {
  method <- match.arg(method, choices = c("first", "last", "average"))
  phylodist <- phylodist(tree_tbl)
  colldist <- colldist(tree_tbl, ...)
  
  mrca <- dplyr::full_join(phylodist, colldist, by = c("tip1", "tip2"))
  
  mrca$mrca <- switch(method,
                      first = (mrca$phylodist-mrca$colldist)/2,
                      last = (mrca$phylodist+mrca$colldist)/2,
                      average = mrca$phylodist/2)
  attr(mrca, "mrca_method") <- method
  return(mrca)
}
