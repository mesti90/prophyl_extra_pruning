#' Are two countries neighbors?
#' 
#' This function uses an external database to decide whether two countries have
#' land borders.
#' 
#' @param x character; country identifier, e.g. name
#' @param y character; country identifier, e.g. name
#' @param from character; identifier type, can be one of \code{"cctld"},
#' \code{"country.name"}, \code{"country.name.de"}, \code{"country.name.fr"},
#' \code{"country.name.it"}, \code{"cowc"}, \code{"cown"}, \code{"dhs"},
#' \code{"ecb"}, \code{"eurostat"}, \code{"fao"}, \code{"fips"}, \code{"gaul"},
#' \code{"genc2c"}, \code{"genc3c"}, \code{"genc3n"}, \code{"gwc"},
#' \code{"gwn"}, \code{"imf"}, \code{"ioc"}, \code{"iso2c"}, \code{"iso3c"},
#' \code{"iso3n"}, \code{"p4c"}, \code{"p4n"}, \code{"un"}, \code{"un_m49"},
#' \code{"unicode.symbol"}, \code{"unhcr"}, \code{"unpd"}, \code{"vdem"},
#' \code{"wb"}, \code{"wb_api2c"}, \code{"wb_api3c"}, \code{"wvs"},
#' \code{"country.name.en.regex"}, \code{"country.name.de.regex"},
#' \code{"country.name.fr.regex"} or \code{"country.name.it.regex"}
#' @references https://github.com/wmgeolab/rgeoboundaries
#' @examples
#' # Hungary and Slovakia have a land border 
#' neighbors("hungary", "slovakia")
#' 
#' # Hungary and Greece do not have a land border
#' neighbors("hungary", "greece")
#' @export
neighbors <- function(x, y, id.type = "country.name") {
  data(country_borders, envir = rlang::current_env())
  iso2c_x <- countrycode::countrycode(
    x, origin = id.type, destination = "iso2c")
  iso2c_y <- countrycode::countrycode(
    y, origin = id.type, destination = "iso2c")
  res <- ifelse(sum(
    country_borders$iso2c1 %in% iso2c_x & country_borders$iso2c2 %in% iso2c_y),
    TRUE,
    FALSE
  )
  return(res)
}
