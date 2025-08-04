#' Get thresholds for non-outlier range based on quantiles
#'
#' This function is a wrapper around the \code{\link{quantile}} function. It can 
#' be used for calculating thresholds below or above which a value should be 
#' considered an outlier. By default, it returns the standard boxplot 
#' thresholds, i.e. the upper threshold is the third quartile plus 1.5 times the 
#' interquartile range (IQR), the lower threshold is the first quartile minus
#' 1.5 times the IQR.
#' @param x numeric; a vector of values
#' @param qlow numeric; lower quantile of a quantile range
#' @param qhigh numeric; higher quantile for a quantile range
#' @param m numeric; range multiplier
#' @param keep character; the value range you would like to keep. Can be
#' \code{"higher"}, \code{"lower"} or \code{"middle"}.
#' @param exclude numeric; values to exclude from the calculating the
#' thresholds.
#' @return a numeric vector with two elements: the lower and the upper 
#' threshold.
#' @examples
#' set.seed(0)
#' data <- rnorm(10000, 0, 1)
#' # Return standard box-whisker outlier thresholds
#' get_non_outlier_range(
#'   data,
#'   qlow = 0.25, 
#'   qhigh = 0.75, 
#'   m = 1.5,
#'   keep = "middle"
#' )
#' # Return lower threshold for values well below the 0.05 quantile
#' get_non_outlier_range(
#'   data,
#'   qlow = 0.05,
#'   qhigh = 1,
#'   m = 1, 
#'   keep = "higher"
#' )
#' # Return upper threshold for values well above the 0.95 quantile
#' get_non_outlier_range(
#'   data,
#'   qlow = 0,
#'   qhigh = 0.95,
#'   m = 1,
#'   keep = "lower"
#' )
#' # Return upper threshold, but exclude values of 0 from calculation
#' get_non_outlier_range(
#'   data,
#'   qlow = 0,
#'   qhigh = 0.95,
#'   m = 1,
#'   keep = "lower",
#'   exclude = 0
#' )
get_non_outlier_range <- function(
    x, 
    qlow = 0.25, 
    qhigh = 0.75, 
    m = 1.5, 
    keep = "middle", 
    exclude = NULL  
  ){
  if (!is.null(exclude)) {
    Qlow <- unname(quantile(x[which(!x %in% exclude)], qlow, na.rm = TRUE))
    Qhigh <- unname(quantile(x[which(!x %in% exclude)], qhigh, na.rm = TRUE))
  } else {
    Qlow <- unname(quantile(x, qlow, na.rm = TRUE))
    Qhigh <- unname(quantile(x, qhigh, na.rm = TRUE))
  }
  lower_threshold <- Qlow - m * (Qhigh-Qlow)
  higher_threshold <- Qhigh + m * (Qhigh-Qlow)
  if (keep == "higher") {
    range <- c(lower_threshold, NA)
  }
  if (keep == "lower") {
    range <- c(NA, higher_threshold)
  }
  if (keep == "middle") {
    range <- c(lower_threshold, higher_threshold)
  }
  return(range)
}
