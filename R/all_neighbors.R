#' List all neighbors of a country
#' 
#' This function uses an external database to list all neighbors of a country.
#' 
#' @param x character; the two character ISO country code of a country 
#' @param borders a data.frame which contains the country borders.
#' @note The database which is used by default does not include e.g. borders of
#' countries that are not universally recognised. However, you can customise
#' your database by adding or removing borders between countries. To do this,
#' use the \code{edit_borders} function and then use the resulting data frame
#' in the \code{borders} argument of this function. 
#' @references https://github.com/wmgeolab/rgeoboundaries
#' @examples
#' # USA has borders with Canada and Mexico
#' all_neighbors("US")
#' 
#' # Kosovo is not recovnised in the default database. As a results, Serbia has
#' # border with Albania
#' all_neighbors("RS")
#' 
#' # However, if we edit the borders database to recognise Kosoveo, Serbia will
#' # no longer have a border with Albania, insted it will have a border with
#' # Kosovo.
#' data(custom_country_borders)
#' updated_borders <- edit_borders(custom_country_borders)
#' all_neighbors("RS", borders = updated_borders)
#' @export
all_neighbors <- function(x, borders = NULL) {
  if (is.null(borders)) {
    data(country_borders, envir = rlang::current_env())
  } else {
    country_borders = borders
  }
  index <- which(country_borders$iso2c1 == x | country_borders$iso2c2 == x)
  iso2c <- sort(unique(unname(unlist(country_borders[index,]))))
  res <- iso2c[-which(iso2c == x)]
  return(res)
}
