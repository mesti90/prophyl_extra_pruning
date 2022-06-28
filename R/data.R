#' Country Borders
#'
#' A data frame of country borders. Each row is a border between two countries.
#'
#' @format A data frame with 728 rows and 2 columns:
#' \describe{
#'   \item{iso2c1}{ISO country code for country 1, two characters}
#'   \item{iso2c2}{ISO country code for country 1, two characters}
#' }
#' @source \url{https://github.com/geodatasource/country-borders}
"country_borders"

#' Custom Country Borders
#'
#' A data frame of custom country borders. Each row is a border between two
#' countries. This data frame includes borders of countries that are not
#' included in the original source e.g. because they are not universally
#' recognised. These rows are then used to modify the original database.
#'
#' @format A data frame with 728 rows and 2 columns:
#' \describe{
#'   \item{iso2c1}{ISO country code for country 1, two characters}
#'   \item{iso2c2}{ISO country code for country 1, two characters}
#'   \item{action}{Whether the row should be added to or removed from database}
#' }
"country_borders"
