#' Compare states between tips
#'
#' This function compares character states for all pairs of tips within a
#' phylogenetic tree.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param state_var character; the name of the variable which stores the
#' character state.
#' @return a tibble with three columns. The third column is 0 if the character
#' states are the same and 1 if they are different.
#' @examples
#' \dontrun {
#' # TODO put a tree_tbl into data so it can be imported.
#' statediff(tree_tbl, state_var = "k_serotype")
#' }
statediff <- function(tree_tbl, state_var) {
  if (state_var %in% names(tree_tbl) == FALSE) {
    stop(state_var, " not found.")
  }
  if (class(tree_tbl[[state_var]]) != "character") {
    stop(state_var, " must be a 'character'.")
  }
  treetips <- tree_tbl[which(grepl("^Node_", tree_tbl$label) == FALSE),]
  statediff <- matrix(1, nrow(treetips), nrow(treetips))
  rownames(statediff) <- treetips$label
  colnames(statediff) <- treetips$label
  for (i in unique(treetips[[state_var]])){
    index = which(treetips[[state_var]] == i)
    statediff[index, index] = 0
  }
  statediff[lower.tri(statediff)] <- NA
  diag(statediff)<-NA
  statediff_df <- dplyr::bind_cols(
    "tip1" = row.names(statediff),
    as.data.frame(statediff)
  )
  statediff_df <- tidyr::pivot_longer(
    statediff_df,
    cols = 2:ncol(statediff_df),
    names_to = "tip2",
    values_to = state_var
  )
  statediff_df <- statediff_df[which(!is.na(statediff_df[[state_var]])),]
  return(statediff_df)
}
