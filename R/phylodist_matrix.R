#' Temporal (cophenetic) distance matrix
#' 
#' This function calculates a matrix which tells the cophenetic distance between
#' pairs of samples. This matrix is one of the "distance" matrices used in risk
#' calculations.
#' @param tree treedater; a dated phylogenetic tree
#' @param df data.frame; a data frame containing sample metadata
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @return a matrix
#' @examples
#' # TODO write an example
#' set.seed(0)
#' tr <- ape::rtree(10)
#' # date tree, define df, calculate matrix
phylodist_matrix <- function(tree,
                             df,
                             focus_on = NULL,
                             focus_by = NULL) {
  # drop tips which cannot be found in assembly table and give a warning
  index <- which(tree$tip.label %in% df$assembly == FALSE)
  if (length(index) > 0) {
    tips_to_drop <- tree$tip.label[index]
    tips_to_drop_collapsed <- paste(tips_to_drop, collapse = ", ")
    tree <- ape::drop.tip(tree, tips_to_drop)
    msg <- paste0(
      "One or more tips could not be found in assembly table and were dropped: ",
      tips_to_drop_collapsed,
      "."
    )
    warning(msg)
  }
  # filter to assemblies that are included in the tree
  index <- which(df$assembly %in% tree$tip.label == FALSE)
  if (length(index) > 0) {
    df <- df[-index, ]
  }
  # calculate cophenetic distance matrix
  phylodist <- ape::cophenetic.phylo(tree)
  # mask matrix if analysis is "focused"
  if (!is.null(focus_on) & !is.null(focus_by)) {
    maskmat <- mask_matrix(df = df, focus_on = focus_on, focus_by = focus_by)
    # reorder maskmat rows and columns to match phylodist rows and columns
    index_mask <- sapply(colnames(phylodist), function(x) {
      which(colnames(maskmat) == x) 
    })
    print(index_mask)
    maskmat_ordered <- maskmat[index_mask, index_mask]
    phylodist <- phylodist * maskmat_ordered
  }
  diag(phylodist) <- NA
  return(phylodist)
}
