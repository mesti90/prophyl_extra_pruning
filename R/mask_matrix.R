#' Mask a distance matrix
#'
#' This helper function is used in distance calculations where we do want to
#' include all tips of a phylogenetic tree and perform pairwise comparisons, 
#' but we also want to define a focus group and exclude comparisons where
#' neither of the two samples belong to the focus group. This function creates
#' a masking matrix that can be used to replace these comparisons with NAs.
#' @param df data.frame; a data frame containing sample metadata
#' @param id_var character; a variable of the df we want to use as row and
#' column names. Must be unique for each row of the data frame.
#' @param focus_by character; a variable of the df we want to focus on
#' @param focus_on character; one or more values of focus_by we want to focus on
#' @return a matrix
#' @examples
#' df <- data.frame(
#'   assembly = LETTERS[1:5],
#'   continent = c("europe", "europe", "asia", "europe", "asia"),
#'   country = c("hungary", "croatia", "china", "serbia", "thailand")
#' )
#' # focus on European samples
#' mask_matrix(
#'   df,
#'   id_var = "assembly",
#'   focus_by = "continent",
#'   focus_on = "europe"
#' )
#' # focus on samples from Croatia or Thailand
#' mask_matrix(
#'   df,
#'   id_var = "assembly",
#'   focus_by = "country",
#'   focus_on = c("croatia", "thailand")
#' )
mask_matrix <- function(df,
                        id_var,
                        focus_by = NULL,
                        focus_on = NULL) {
  
  # validate id_var
  id_var <- validate_id_var(df, id_var)
  # validate focus_by and focus_on
  focus <- validate_focus(df, focus_by, focus_on)
  focus_by <- focus$focus_by
  focus_on <- focus$focus_on
  
  mask <- matrix(1, nrow(df), nrow(df))
  rownames(mask) <- df[[id_var]]
  colnames(mask) <- df[[id_var]]
  
  if (!is.null(focus_by) & !is.null(focus_on)) {
    index <- which(df[[focus_by]] %in% focus_on == FALSE)
    if (length(index) == 0) {
      warning("Focus group contains all isolates.")
    } else {
      mask[index, index] <- NA
    }
  }
  
  return(mask)
}
