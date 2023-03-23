#' Neighbors matrix
#' 
#' This function calculates a matrix which tells whether pairs of samples were
#' collected from neighboring countries. This matrix is one of the "distance"
#' matrices used in risk calculations.
#' @param df data.frame; a data frame containing sample metadata
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @return a matrix
#' @note The function requires ISO 3166-1 alpha-2 codes for countries. These 
#' can be found using \code{countrycode::countrycode()}
#' @examples
#' # country code for china
#' countrycode::countrycode(
#'   "china", origin = "country.name", destination = "iso2c")
#' df <- data.frame(
#'   assembly = LETTERS[1:5],
#'   continent = c("europe", "europe", "asia", "europe", "asia"),
#'   country = c("hungary", "germany", "china", "serbia", "laos"),
#'   country_iso2c = c("HU", "DE", "CN", "RS", "LA")
#' )
#' # no focus
#' neighbors_matrix(df)
#' # focus on Asian samples
#' neighbors_matrix(df, focus_by = "continent", focus_on = "asia")
neighbors_matrix <- function(df,
                             focus_on = NULL,
                             focus_by = NULL) {
  if ("country_iso2c" %in% names(df) == FALSE) {
    stop("Required variable 'country_iso2c' is missing from df.")
  }
  if (!is.null(focus_by)) {
    focus_by <- match.arg(focus_by, choices = names(df))
  }
  if (!is.null(focus_on)) {
    focus_on <- match.arg(focus_on, choices = unique(df[[focus_by]]))
  }
  neighbors <- matrix(0, nrow(df), nrow(df))
  for (i in unique(df$country_iso2c)){
    index1 <- which(df$country_iso2c == i)
    # apply custom country borders
    data("custom_country_borders")
    borders <- edit_borders(custom_country_borders)
    index2 <- which(df$country_iso2c %in% all_neighbors(i, borders = borders))
    
    neighbors[index1, index2] <- 1
    neighbors[index2, index1] <- 1
  }
  if (!is.null(focus_on) & !is.null(focus_by)) {
    maskmat <- mask_matrix(df = df, focus_on = focus_on, focus_by = focus_by)
    neighbors <- neighbors * maskmat
  }
  diag(neighbors)<-NA
  rownames(neighbors) <- df$assembly
  colnames(neighbors) <- df$assembly
  return(neighbors)
}
