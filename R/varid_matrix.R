#' Variable identity matrix
#' 
#' This function calculates a matrix which tells whether a selected variable is 
#' the same for pairs of samples, e.g. whether pairs fo samples were collected
#' in the same country. This matrix is one of the "distance" matrices used in
#' risk calculations.
#' @param df data.frame; a data frame containing sample metadata
#' @param id_var character; a variable of the df we want to use as row and
#' column names. Must be unique for each row of the data frame.
#' @param var character; a variable of the df to assess
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @return a matrix
#' @examples
#' df <- data.frame(
#'   assembly = LETTERS[1:5],
#'   continent = c("europe", "europe", "asia", "europe", "asia"),
#'   country = c("hungary", "hungary", "china", "serbia", "china"),
#'   city = c("budapest", "szeged", "beijing", "belgrade", "beijing")
#' )
#' # which samples come from the same country
#' varid_matrix(df, id_var = "assembly", var = "country")
#' # same but focus on European samples
#' varid_matrix(
#'   df, 
#'   id_var = "assembly", 
#'   var = "country", 
#'   focus_by = "continent", 
#'   focus_on = "europe"
#'  )
#' # which samples have the same mlst
#' varid_matrix(df, id_var = "assembly", var = "city")
varid_matrix <- function(df,
                         id_var,
                         var = "country",
                         focus_on = NULL,
                         focus_by = NULL) {
  
  # validate id_var
  id_var <- validate_id_var(df, id_var)
  # validate var
  var <- match.arg(var, choices = names(df))
  # validate focus_by and focus_on
  focus <- validate_focus(df, focus_by, focus_on)
  focus_by <- focus$focus_by
  focus_on <- focus$focus_on
  varid <- matrix(0, nrow(df), nrow(df))
  rownames(varid) <- df[[id_var]]
  colnames(varid) <- df[[id_var]]
  for (i in unique(df[[var]])){
    if (is.na(i)) {
      index <- which(is.na(df[[var]]))
      varid[index, ] <- NA
      varid[, index] <- NA
    } else {
      index <- which(df[[var]] == i)
      varid[index, index] <- 1
    }
  }
  maskmat <- mask_matrix(
    df = df,
    id_var = id_var,
    focus_on = focus_on,
    focus_by = focus_by
  )

  if (any(row.names(maskmat) != row.names(varid))) {
    stop("Matrix row and/or columns are mixed up.")
  }
  if (any(colnames(maskmat) != colnames(varid))) {
    stop("Matrix row and/or columns are mixed up.")
  }

  varid <- varid * maskmat
  diag(varid)<-NA
  return(varid)
}
