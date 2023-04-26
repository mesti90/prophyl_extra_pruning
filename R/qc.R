#' Check whether a set of values fit into a quality range
#' 
#' @param x numeric; a vector of values
#' @param qlow numeric; lower quantile of a quantile range
#' @param qhigh numeric; higher quantile for a quantile range
#' @param m numeric; range multiplier
#' @param keep character; the value range you would like to keep. Can be
#' \code{"higher"}, \code{"lower"} or \code{"middle"}.
#' @param exclude_zero logical; 
#' @examples
#' set.seed(0)
#' data <- rnorm(100, 0, 1)
#' # Calculate interquartile range (IQR), keep middle range
#' qc(data, qlow = 0.25, qhigh = 0.75, m = 1.5, keep = "middle")
#' # Keep high values (flag values well below the 0.05 quantile)
#' qc(data, qlow = 0.05, qhigh = 1, m = 1, keep = "higher")
#' # Keep low value (flag values well above the 0.95 quantile)
#' qc(data, qlow = 0, qhigh = 0.95, m = 1, keep = "lower")
#' @export
qc <- function(x, qlow = 0, qhigh = 1, m, keep, exclude) {
  Qlow <- unname(quantile(x, qlow, na.rm = TRUE))
  Qhigh <- unname(quantile(x, qhigh, na.rm = TRUE))
  lower_threshold <- Qlow - m * (Qhigh-Qlow)
  higher_threshold <- Qhigh + m * (Qhigh-Qlow)
  if (keep == "higher") {
    out <- ifelse(
      !is.na(x),
      x >= lower_threshold,
      NA
    )
  }
  if (keep == "lower") {
    
    out <- ifelse(
      !is.na(x),
      x < higher_threshold,
      NA
    )
  }
  if (keep == "middle") {
    out <- ifelse(
      !is.na(x), 
      x >= lower_threshold & x < higher_threshold,
      NA
    )
  }
  out <- list(
    value = out,
    thresholds = c(lower_threshold, higher_threshold)
  )
  return(out)
}
# NOTE in this application it may be easier to just calculate Q3 + 3 * (Q3 - Q1)
# instead of using a separate function