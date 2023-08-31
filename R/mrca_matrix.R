#' Temporal distance (MRCA) matrix
#' 
#' This function calculates a matrix which tells the date of the most recent
#' common ancestor (MRCA) between pairs of samples. This matrix is one of the
#' "distance" matrices used in risk calculations.
#' @param phylodist phylodist; a cophenetic distance matrix
#' @param colldist colldist; a collection date distance matrix
#' @param force_nonnegative logical; should negative values be converted to 0?
#' @return a matrix
#' @seealso [phylodist_matrix()], [tempdist_matrix()] 
#' @examples
#' df <- data.frame(
#'   assembly = LETTERS[1:5],
#'   continent = c("europe", "europe", "asia", "europe", "asia"),
#'   country = c("hungary", "germany", "china", "hungary", "laos"),
#'   collection_date = c(
#'     "2000-01-01", "2001-01-01", "2005-01-01", "2000-06-15", "2005-06-15" 
#'   )
#' )
#' colldist <- tempdist_matrix(
#'   df, id_var = "assembly", date_var = "collection_date")
#' set.seed(0)
#' tr <- ape::rtree(5, tip.label = df$assembly)
#' phylodist <- phylodist_matrix(tree = tr, df = df,id_var = "assembly")
#' mrca_matrix(phylodist, colldist)
mrca_matrix <- function(phylodist,
                        colldist,
                        force_nonnegative = TRUE) {
  if (any(row.names(phylodist) != row.names(colldist))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  if (any(colnames(phylodist) != colnames(colldist))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  # definition:
  # mrca is the distance of the common ancestor from the older sample
  mrca <- (phylodist - colldist)/2
  # alternative definition:
  # mrca is the distance of the common ancestor from the more recent sample
  # TODO discuss because the alternative makes more sense also we might be able
  # to eliminated the sampling window with it?
  # mrca <- (phylodist + colldist)/2
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
