#' Shrink matrix
#' 
#' This function shrinks a matrix to have the same rows and columns as a smaller
#' matrix and also in the same order which is important for later calculations.
#' @param matrix matrix; a matrix to shrink
#' @param small_matrix matrix; a matrix to shrink to
#' @return a matrix
shrink_matrix <- function(matrix, small_matrix) {
  index <- unname(sapply(colnames(small_matrix), function(x) {
    which(colnames(matrix) == x)
  }))
  shrinked_matrix <- matrix[index, index]
  if (any(row.names(shrinked_matrix) != row.names(small_matrix))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  if (any(colnames(shrinked_matrix) != colnames(small_matrix))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  return(shrinked_matrix)
}
