#' Temporal distance (MRCA) matrix
#' 
#' This function calculates a matrix which tells the date of the most recent
#' common ancestor (MRCA) between pairs of samples. This matrix is one of the
#' "distance" matrices used in risk calculations.
#' @param phylodist phylodist; a cophenetic distance matrix
#' @param colldist colldist; a collection date distance matrix
#' @param force_nonnegative logical; should negative values be converted to 0?
#' @return a matrix
#' @seealso [phylodist_matrix()], [colldist_matrix()] 
#' @examples
#' # TODO
mrca_matrix <- function(phylodist,
                        colldist,
                        force_nonnegative = TRUE) {
  if (all(colnames(phylodist) %in% colnames(colldist)) == FALSE) {
    stop()
  }
  if (all(colnames(colldist) %in% colnames(phylodist)) == FALSE) {
    stop()
  }
  # reorder colldist rows and columns to match phylodist rows and columns
  index_colldist <- unname(sapply(colnames(phylodist), function(x) {
    which(colnames(colldist) == x)
  }))
  colldist_ordered <- colldist[index_colldist, index_colldist]
  # applied mrca definition: distance from the older sample 
  mrca <- (phylodist - colldist_ordered)/2
  if (force_nonnegative == TRUE) {
    index <- which(mrca < 0)
    if (length(index) > 0) {
      mrca[which(mrca < 0)] <- 0
      msg <- paste0(
        "MRCA value for ", length(index)/2, " pairs was below 0. ",
        "These values were set to 0."
      )
      warning(msg)
    }
  }
  return(mrca)
}
