#' Temporal (sampling) distance matrix
#' 
#' This function calculates a matrix which tells the time difference between
#' sampling dates for pairs of samples. This matrix is one of the "distance"
#' matrices used in risk calculations.
#' @param df data.frame; a data frame containing sample metadata
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @param estimate_dates character; strategy for estimating uncertain dates.
#' Either \code{"asis"} (no estimation), \code{"upper"}, \code{"middle"}, or 
#' \code{"lower"}.
#' @seealso [date_upper()],[date_middle()],[date_lower()]
#' @return a matrix
#' @examples
#' df <- data.frame(
#'   assembly = LETTERS[1:5],
#'   continent = c("europe", "europe", "asia", "europe", "asia"),
#'   country = c("hungary", "germany", "china", "hungary", "laos"),
#'   collection_date = c(
#'     "2000-01-01", "2001-01-01", "2005-01-01", "2000-06-15", "2005-06-15" 
#'   )
#' )
#' # no focus
#' colldist_matrix(df)
#' # focus on European samples
#' colldist_matrix(df, focus_by = "continent", focus_on = "asia")
colldist_matrix <- function(df,
                            focus_on = NULL,
                            focus_by = NULL,
                            estimate_dates = "asis") {
  if ("collection_date" %in% names(df) == FALSE) {
    "Required variable 'collection_date' is missing from df."
  }
  estimate_dates <- match.arg(estimate_dates, choices = c(
    "asis", "upper", "middle", "lower" 
  ))
  dates <- switch(
    estimate_dates,
    asis = as.Date(df$collection_date),
    upper = date_upper(df$collection_date),
    middle = date_middle(df$collection_date),
    lower = date_lower(df$collection_date)
  )
  dates <- unname(lubridate::decimal_date(dates))
  colldist <- round(abs(outer(dates, dates, "-")),2)
  if (!is.null(focus_on) & !is.null(focus_by)) {
    maskmat <- mask_matrix(df = df, focus_on = focus_on, focus_by = focus_by)
    colldist <- colldist * maskmat
  }
  diag(colldist) <- NA
  rownames(colldist) <- df$assembly
  colnames(colldist) <- df$assembly
  return(colldist)
}
