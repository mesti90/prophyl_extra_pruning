#' Distances between sample collection dates
#' 
#' This functions calculates distances between sample collection dates for all
#' pairs of tips within a phylogenetic tree.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param date_var character; the name of the variable which stores sample
#' collection dates.
#' @return a tibble with three columns.
#' @examples
#' \dontrun{
#' # TODO put a tree_tbl into data so it can be imported.
#' colldist(tree_tbl)
#' }
colldist <- function(tree_tbl, date_var = "collection_day") {
  #TODO: validate tree_tbl input
  if (date_var %in% names(tree_tbl) == FALSE) {
    stop(date_var, " not found.")
  }
  if (class(tree_tbl[[date_var]]) != "Date") {
    stop(date_var, " must be a 'Date'.")
  }
  treetips <- tree_tbl[which(grepl("^Node_", tree_tbl$label) == FALSE),]
  dates <- unname(lubridate::decimal_date(treetips[[date_var]]))
  colldist = round(abs(outer(dates, dates, "-")),2)
  colldist[lower.tri(colldist)] <- NA
  diag(colldist) <- NA
  rownames(colldist) <- treetips$label
  colnames(colldist) <- treetips$label
  colldist_df <- dplyr::bind_cols(
    "tip1" = row.names(colldist),
    as.data.frame(colldist)
  )
  colldist_df <- tidyr::pivot_longer(
    colldist_df,
    cols = 2:ncol(colldist_df),
    names_to = "tip2",
    values_to = "colldist"
  )
  colldist_df <- colldist_df[which(!is.na(colldist_df$colldist)),]
  return(colldist_df)
}
