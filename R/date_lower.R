#' Convert date to numeric format for phylogenetic dating
#'
#' This is a utility function that converts a date to a numeric value. This
#' format is commonly used in phylogenetic dating. If date is incomplete the
#' function returns the lower end of the interval.
#' @param dates character; a full or partial date in "YYYY-MM_DD" like format.
#' @param out_format character; output format, either \code{"decimal"} or
#' \code{"date"}.
#' @examples
#' date_lower("2021-08-11")
#' date_lower("1988-03")
#' date_lower("2000")
#' date_lower("1900/1948")
#' date_lower("2003-03/2005-02")
#' @export
date_lower <- function(dates, out_format = "decimal") {
  
  out_format <- match.arg(out_format, choices = c("decimal", "date"))
  
  foo <- function(x) {
  if (is.na(x)) return(NA)
  year = suppressWarnings(
    as.numeric(stringi::stri_sub(x, from = 1, to = 4))
  )
  if (is.na(year)) {
    msg <- paste0("Date conversion failed: '", x, "'. Returning NA.")
    warning(msg)
    return(NA)
  }
  date_elements <- strsplit(x, split = "-")[[1]]
  year <- date_elements[1]
  year_start <- as.Date(paste0(year, "-01-01"))
  year_end <- as.Date(paste0(year, "-12-31"))
  if (length(date_elements) == 3) {
    date <- as.Date(x, format = "%Y-%m-%d")
  }
  if (length(date_elements) == 2) {
    month <- date_elements[2]
    date <- as.Date(paste(year, month, "01", sep = "-"))
  }
  if (length(date_elements) == 1) {
    date <- as.Date(paste(year, "01", "01", sep = "-"))
  }
  date_out <- as.numeric(date-year_start)/as.numeric(year_end-year_start+1)
  # 1st january may be problematic
  date_out <- round(1000*date_out, 0)
  if(nchar(date_out) == 2){
    date_out <- paste0("0", date_out)
  }
  if(nchar(date_out) == 1) {
    date_out <- paste0("00", date_out)
  }
  date_out <- paste0(year, ".", date_out)

  return(date_out)
  }
  
  bar <- function(x) {
    candidates <- strsplit(x, split = "/")[[1]]
    lows <- sapply(candidates, foo)
    return(min(lows))
  }
  
out <- unname(sapply(dates, function(x) try(bar(x), silent = TRUE)))
out <- as.numeric(out)

if (out_format == "date") {
  out <- as.Date(lubridate::date_decimal(out), format = "%Y-%m-%d")
}

return(out)
}
