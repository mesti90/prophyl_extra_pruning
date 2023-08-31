#' Simulate random dates
#' 
#' This function returns exact dates and simulates random dates for incomplete
#' dates. For incomplete dates, the function will attempt to guess the lower and
#' upper boundaries from the input string.
#' @param dates character; a vector of dates. Can be exact, e.g. 2000-01-01, or
#' incomplete, e.g. 2000-01 or 2000.
#' @param out_format character; output format, either \code{"decimal"} or
#' \code{"date"}.
#' @return a vector of dates either in decimal format. If the input is complete,
#' the output will be the same as the input. If the input is incomplete, the
#' output will be a random date between the lower and upper boundaries of the
#' input.
#' @examples
#' date_runif("2000-01-01")
#' date_runif("2000-01")
#' date_runif("2000")
#' @export
date_runif <- function(dates, out_format = "decimal") {
  
  out_format <- match.arg(out_format, choices = c("decimal", "date"))
  
  dates <- as.character(dates)
  foo <- function(x) {
    if (is.na(x)) return(NA)
    lower <- date_lower(x)
    upper <- date_upper(x)
    if (!is.na(lower) & !is.na(upper)) {
      runif(1, min = lower, max = upper)
    } else {
      msg <- paste0("Cannot simulate date for ", x, ". Returning NA.")
      warning(msg)
      return(NA)
    }
  }
  out <- unname(sapply(dates, foo))
  
  if (out_format == "date") {
    out <- as.Date(lubridate::date_decimal(out), format = "%Y-%m-%d")
  }
  
  return(out)
}
