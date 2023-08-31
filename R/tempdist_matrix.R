#' Temporal distance matrix
#' 
#' This function calculates a matrix which tells the time difference between
#' sampling dates for pairs of samples. This matrix is one of the "distance"
#' matrices used in risk calculations.
#' @param df data.frame; a data frame containing sample metadata
#' @param id_var character; a variable of the df we want to use as row and
#' column names. Must be unique for each row of the data frame.
#' @param date_var character; a variable of the df we want to assess which
#' stores dates.
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @param estimate_dates character; strategy for estimating uncertain dates.
#' Either \code{"asis"} (no estimation), \code{"upper"}, \code{"middle"}, 
#' \code{"lower"} or \code{"runif"}.
#' @seealso [date_upper()],[date_middle()],[date_lower()], [date_runif()]
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
#' tempdist_matrix(df, id_var = "assembly", date_var = "collection_date")
#' # focus on European samples
#' tempdist_matrix(
#'   df,
#'   id_var = "assembly",
#'   date_var = "collection_date",
#'   focus_by = "continent",
#'   focus_on = "asia"
#' )
tempdist_matrix <- function(df,
                            id_var,
                            date_var, 
                            focus_on = NULL,
                            focus_by = NULL,
                            estimate_dates = "asis") {
  # validate id_var
  id_var <- validate_id_var(df, id_var)
  # validate date_var
  date_var <- match.arg(date_var, choices = names(df))
  # validate focus_by and focus_on
  focus <- validate_focus(df, focus_by, focus_on)
  focus_by <- focus$focus_by
  focus_on <- focus$focus_on
  # validate estimate_dates
  estimate_dates <- match.arg(estimate_dates, choices = c(
    "asis", "upper", "middle", "lower", "runif" 
  ))

  dates <- switch(
    estimate_dates,
    asis = as.Date(df[[date_var]]),
    upper = date_upper(df[[date_var]], out_format = "date"),
    middle = date_middle(df[[date_var]], out_format = "date"),
    lower = date_lower(df[[date_var]], out_format = "date"),
    runif = date_runif(df[[date_var]], out_format = "date")
  )
  dates <- unname(lubridate::decimal_date(dates))

  colldist <- round(abs(outer(dates, dates, "-")),2)
  rownames(colldist) <- df[[id_var]]
  colnames(colldist) <- df[[id_var]]

  maskmat <- mask_matrix(
      df = df,
      id_var = id_var,
      focus_on = focus_on,
      focus_by = focus_by
  )
   
  if (any(row.names(maskmat) != row.names(colldist))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  if (any(colnames(maskmat) != colnames(colldist))) {
    stop("Matrix row and/or columns are mixed up.")
  }

  colldist <- colldist * maskmat
  
  diag(colldist) <- NA

  return(colldist)
}
