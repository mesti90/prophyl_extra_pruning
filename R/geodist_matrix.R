#' Geographic distance matrix
#' 
#' This function calculates a matrix which tells the geographic distance in km
#' between pairs of samples. This matrix is one of the "distance" matrices used
#' in risk calculations.
#' @param df data.frame; a data frame containing sample metadata
#' @param id_var character; a variable of the df we want to use as row and
#' column names. Must be unique for each row of the data frame.
#' @param lat_var character; a variable of the df we want to assess which
#' stores latitude coordinate.
#' @param lon_var character; a variable of the df we want to assess which
#' stores longitude coordinate.
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
#' geodist_matrix(df, id_var = "assembly")
#' # focus on Asian samples
#' geodist_matrix(
#'   df,
#'   id_var = "assembly",
#'   focus_by = "continent",
#'   focus_on = "asia"
#' )
geodist_matrix <- function(df,
                           id_var,
                           lat_var = "lat",
                           lon_var = "lon",
                           focus_on = NULL,
                           focus_by = NULL) {

  # validate id_var
  id_var <- validate_id_var(df, id_var)
  # validate lat_var
  lat_var <- match.arg(lat_var, choices = names(df))
  # validate lon_var
  lon_var <- match.arg(lon_var, choices = names(df))
  # validate focus_by and focus_on
  focus <- validate_focus(df, focus_by, focus_on)
  focus_by <- focus$focus_by
  focus_on <- focus$focus_on

  geodist <- matrix(NA, nrow(df), nrow(df))
  rownames(geodist) <- df[[id_var]]
  colnames(geodist) <- df[[id_var]]
  
  indices <- 1:nrow(df)
  for (i in 1:nrow(df)) {
    lat1 <- df[[lat_var]][i]
    lon1 <- df[[lon_var]][i]
    # if coordinates are identical set distance to 0.
    index <- which(df[[lat_var]] == lat1 & df[[lon_var]] == lon1)
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
        new <- which(df[[lon_var]] == s$lon[j] & df[[lat_var]] == s$lat[j])
        geodist[index, new] <- geodist_km
        geodist[new, index] <- geodist_km
      }
    }
  }
  
  maskmat <- mask_matrix(
    df = df,
    id_var = id_var,
    focus_on = focus_on,
    focus_by = focus_by
  )

  if (any(row.names(maskmat) != row.names(geodist))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  if (any(colnames(maskmat) != colnames(geodist))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  
  geodist <- geodist * maskmat
  
  diag(geodist) <- NA

  return(geodist)
}
