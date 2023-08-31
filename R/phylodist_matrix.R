#' Temporal (cophenetic) distance matrix
#' 
#' This function calculates a matrix which tells the cophenetic distance between
#' pairs of samples. This matrix is one of the "distance" matrices used in risk
#' calculations.
#' @param tree treedater; a dated phylogenetic tree
#' @param df data.frame; a data frame containing sample metadata
#' @param id_var character; a variable of the df we want to use as row and
#' column names. Must be unique for each row of the data frame. Must match
#' tree tip labels.
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @return a matrix
#' @examples
#' df <- data.frame(
#'   assembly = LETTERS[1:5],
#'   continent = c("europe", "europe", "asia", "europe", "asia"),
#'   country = c("hungary", "hungary", "china", "serbia", "china"),
#'   mlst = c("ST1", "ST2", "ST1", "ST2", "ST3")
#' )
#' set.seed(0)
#' tr <- ape::rtree(6, tip.label = df$assembly)
#' phylodist_matrix(tree = tr, df = df,id_var = "assembly")
phylodist_matrix <- function(tree,
                             df,
                             id_var, 
                             focus_on = NULL,
                             focus_by = NULL) {
  
  # validate tree
  if (length(tree$tip.label) > length(unique(tree$tip.label))) {
    stop("Tree tip labels must be unique.")
  }
  # validate id_var
  id_var <- validate_id_var(df, id_var)
  # validate focus_by and focus_on
  focus <- validate_focus(df, focus_by, focus_on)
  focus_by <- focus$focus_by
  focus_on <- focus$focus_on
  
  # drop tips which cannot be found in assembly table and give a warning
  index <- which(tree$tip.label %in% df[[id_var]] == FALSE)
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
  index <- which(df[[id_var]] %in% tree$tip.label == FALSE)
  if (length(index) > 0) {
    df <- df[-index, ]
  }
  # calculate cophenetic distance matrix
  phylodist <- ape::cophenetic.phylo(tree)
  # calculate masking matri
  maskmat <- mask_matrix(
    df = df,
    id_var = id_var,
    focus_on = focus_on,
    focus_by = focus_by
  )
  # reorder phylodist to match mask_matrix
  index <- sapply(colnames(maskmat), function(x) {
    which(colnames(phylodist) == x) 
  })
  phylodist <- phylodist[index,index]
  # test whether rownames and colnames match between the two matrices
  if (any(row.names(maskmat) != row.names(phylodist))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  if (any(colnames(maskmat) != colnames(phylodist))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  phylodist <- phylodist * maskmat
  diag(phylodist) <- NA
  return(phylodist)
}
