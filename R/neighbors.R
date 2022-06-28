#' Are two countries neighbors?
#' 
#' This function uses an external database to decide whether two countries have
#' land borders.
#' 
#' @param x character; a vector of two characters ISO country codes 
#' @param y character; a vector of two characters ISO country codes
#' @param custom_borders an object of class \code{custom_borders} which
#' modifies the default database.
#' @note The database which is used by default does not include e.g. borders of
#' countries that are not universally recognised. However, you can customise
#' your database by adding or removing borders between countries. The
#' \code{custom_borders} argument lets you specify such database.
#' @references https://github.com/wmgeolab/rgeoboundaries
#' @examples
#' # Hungary and Slovakia have a land border 
#' neighbors("HU", "SK")
#' 
#' # Hungary and Greece do not have a land border
#' neighbors("HU", "GR")
#' @export
neighbors <- function(x, y, custom_borders = NULL) {
  # TODO: too slow when many comparions required -> speed up
  data(country_borders, envir = rlang::current_env())
  if (!is.null(custom_borders)) {
    country_borders <- edit_borders(custom_borders)
  }
  res <- unname(mapply(function(a, b) {
   res <-  sum(country_borders$iso2c1 %in% a & country_borders$iso2c2 %in% b)
   res <- ifelse(res > 0, TRUE, FALSE)
   }, x, y))
  return(res)
}
