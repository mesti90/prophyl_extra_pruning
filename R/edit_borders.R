edit_borders <- function(custom_borders) {
  data(country_borders, envir = rlang::current_env())
  # TODO: add validation for custom borders here
  index_add <- which(custom_country_borders$action == "add")
  if (length(index_add) > 0) {
    country_borders <- dplyr::bind_rows(
      country_borders, custom_country_borders[index_add, c(1,2)]
    )
  }
  index_remove <- which(custom_country_borders$action == "remove")
  if (length(index_remove) > 0) {
    for (i in index_remove) {
      idx <- which(
        country_borders$iso2c1 == custom_country_borders$iso2c1[i] &
          country_borders$iso2c2 == custom_country_borders$iso2c2[i])
      if (length(idx) > 0) {
        country_borders <- country_borders[-idx,]
      }
    }
  }
  return(country_borders)
}