#' Neighbors matrix
#' 
#' This function calculates a matrix which tells whether pairs of samples were
#' collected from neighboring countries. This matrix is one of the "distance"
#' matrices used in risk calculations.
#' @param df data.frame; a data frame containing sample metadata
#' @param id_var character; a variable of the df we want to use as row and
#' column names. Must be unique for each row of the data frame.
#' @param iso2c_var character; a variable of the df that contains ISO 3166-1
#' alpha-2 codes for countries. These can be found using e.g.
#' \code{countrycode::countrycode()}.
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @return a matrix
#' @examples
#' # country code for china
#' countrycode::countrycode(
#'   "china", origin = "country.name", destination = "iso2c")
#' df <- data.frame(
#'   assembly = LETTERS[1:7],
#'   continent = c(
#'     "europe", "europe", "asia", "europe", "asia", "europe", "europe"),
#'   country = c(
#'     "hungary", "germany", "china", "serbia", "laos", "albania", "kosovo"),
#'   country_iso2c = c("HU", "DE", "CN", "RS", "LA", "AL", "XK")
#' )
#' # no focus, use default borders
#' neighbors_matrix(df, id_var = "assembly")
#' # no focus, custom borders
#' # load data set with requested changes to the default borders
#' data(custom_country_borders)
#' # update borders and calculate neighbors matrix 
#' updated_borders <- edit_borders(custom_country_borders)
#' neighbors_matrix(df, id_var = "assembly", country_borders = updated_borders)
#' # focus on Asian samples
#' neighbors_matrix(
#'   df,
#'   id_var = "assembly",
#'   focus_by = "continent", 
#'   focus_on = "asia"
#' )
neighbors_matrix <- function(df,
                             id_var,
                             iso2c_var = "country_iso2c",
                             country_borders = "default",
                             focus_on = NULL,
                             focus_by = NULL) {
  
  # validate id_var
  id_var <- validate_id_var(df, id_var)
  # validate iso2c_var
  iso2c_var <- match.arg(iso2c_var, choices = names(df))
  # check whether these are all in database??
  # validate country borders
  if (length(country_borders) == 1 && country_borders == "default") {
    data(country_borders, envir = rlang::current_env())
  }
  all_iso2c <- unique(c(country_borders$iso2c1, country_borders$iso2c2))
  index <- which(!df[[iso2c_var]] %in% all_iso2c)
  if (length(index) > 0) {
    missing_iso2c_collapsed <- paste(
      sort(unique(df[[iso2c_var]][index])), collapse = ", ")
    msg = paste(
    "Some country codes could not be found in the border database: ",
    missing_iso2c_collapsed
    )
    warning(msg)
  }
  # validate focus_by and focus_on
  focus <- validate_focus(df, focus_by, focus_on)
  focus_by <- focus$focus_by
  focus_on <- focus$focus_on

  neighbors <- matrix(0, nrow(df), nrow(df))
  rownames(neighbors) <- df[[id_var]]
  colnames(neighbors) <- df[[id_var]]
  
  borders <- country_borders
  
  for (i in unique(df[[iso2c_var]])){
    if (is.na(i)) {
      index1 <- which(is.na(df[[iso2c_var]]))
      neighbors[index1, ] <- NA
      neighbors[, index1] <- NA
    } else {
      index1 <- which(df[[iso2c_var]] == i)
      if (!i %in% all_iso2c) {
        neighbors[index1, ] <- NA
        neighbors[, index1] <- NA
      } else {
        index2 <- which(df[[iso2c_var]] %in% all_neighbors(i, borders = borders))
        if (length(index2) > 0) {
          neighbors[index1, index2] <- 1
          neighbors[index2, index1] <- 1
        }
      }
    }
  }
  
  maskmat <- mask_matrix(
    df = df,
    id_var = id_var,
    focus_on = focus_on,
    focus_by = focus_by
  )
  
if (any(row.names(maskmat) != row.names(neighbors))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  if (any(colnames(maskmat) != colnames(neighbors))) {
    stop("Matrix row and/or columns are mixed up.")
  }

  neighbors <- neighbors * maskmat
  
  diag(neighbors)<-NA

  return(neighbors)
}
