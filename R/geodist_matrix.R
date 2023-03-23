#' Geographic distance matrix
#' 
#' This function calculates a matrix which tells the geographic distance in km
#' between pairs of samples. This matrix is one of the "distance" matrices used
#' in risk calculations.
#' @param df data.frame; a data frame containing sample metadata
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @return a matrix
#' @examples
#' df <- data.frame(
#'   assembly = LETTERS[1:5],
#'   continent = c("europe", "europe", "asia", "europe", "asia"),
#'   country = c("hungary", "germany", "china", "hungary", "laos"),
#'   lat = c(47.16249, 51.16569, 35.86166, 47.16249, 19.85627),
#'   lon = c(19.50330, 10.451526, 104.19540, 19.50330, 102.4955)
#' )
#' # no focus
#' geodist_matrix(df)
#' # focus on Asian samples
#' geodist_matrix(df, focus_by = "continent", focus_on = "asia")
geodist_matrix <- function(df,
                           focus_on = NULL,
                           focus_by = NULL) {
  if ("lat" %in% names(df) == FALSE) {
    stop("Required variable 'lat' is missing from df.")
  }
  if ("lon" %in% names(df) == FALSE) {
    stop("Required variable 'lon' is missing from df.")
  }
  if (!is.null(focus_by)) {
    focus_by <- match.arg(focus_by, choices = names(df))
  }
  if (!is.null(focus_on)) {
    focus_on <- match.arg(focus_on, choices = unique(df[[focus_by]]))
  }
  geodist <- matrix(NA, nrow(df), nrow(df))
  indices <- 1:nrow(df)
  for (i in 1:nrow(df)) {
    lat1 <- df$lat[i]
    lon1 <- df$lon[i]
    # if coordinates are identical set distance to 0.
    index <- which(df$lat == lat1 & df$lon == lon1)
    geodist[i, index] = 0
    geodist[index, i] = 0
    index_test <- which(is.na(geodist[i, ]))
    if (length(index_test) > 0) {
      s <- df[index_test, which(names(df) %in% c("lat", "lon"))]
      s <- dplyr::distinct(s)
      for (j in 1:nrow(s)) {
        geodist_km <- round(geosphere::distHaversine(
          p1 = c(lon1, lat1),
          p2 = c(s$lon[j], s$lat[j])
        )/1000, 0)
        new <- which(df$lon == s$lon[j] & df$lat == s$lat[j])
        geodist[index, new] <- geodist_km
        geodist[new, index] <- geodist_km
      }
    }
  }
  if (!is.null(focus_on) & !is.null(focus_by)) {
    maskmat <- mask_matrix(df = df, focus_on = focus_on, focus_by = focus_by)
    geodist <- geodist * maskmat
  }
  diag(geodist) <- NA
  rownames(geodist) <- df$assembly
  colnames(geodist) <- df$assembly
  return(geodist)
}
