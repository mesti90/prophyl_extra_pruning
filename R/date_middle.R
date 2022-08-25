#' Convert date to numeric format for phylogenetic dating
#'
#' This is a utility function that converts a date to a numeric value. This
#' format is commonly used in phylogenetic dating. If date is incomplete the
#' function returns the middle of the interval.
#' @param dates character; a full or partial date in "YYYY-MM_DD" like format.
#' @examples
#' date_middle("2021-18-11")
#' date_middle("1988-03")
#' date_middle("2000")
#' @export
date_middle <- function(dates){
  foo <- function(x){
    if(is.na(x)) return(NA)
    year = suppressWarnings(
      as.numeric(stringi::stri_sub(x, from = 1, to = 4))
    )
    if (is.na(year)) {
      return(NA)
    }
    date_elements <- strsplit(x, split = "-")[[1]]
    year <- date_elements[1]
    if (length(date_elements) == 1){
      date <- paste(year, "06","15", sep = "-")
    }
    if (length(date_elements) == 2){
      month <- date_elements[2]
      date <- paste(year, month, "15", sep = "-")
    }
    if (length(date_elements) == 3){
      date <- x
    }
    return(date)
  }
  out <- as.Date(sapply(dates, foo))
  return(out)
}
