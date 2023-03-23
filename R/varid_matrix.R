#' Variable identity matrix
#' 
#' This function calculates a matrix which tells whether a selected variable is 
#' the same for pairs of samples, e.g. whether pairs fo samples were collected
#' in the same country. This matrix is one of the "distance" matrices used in
#' risk calculations.
#' @param df data.frame; a data frame containing sample metadata
#' @param var character; a variable of the df to assess
#' @param focus_by character; a variable of the df to focus on
#' @param focus_on character; one or more values of focus_by to focus on
#' @return a matrix
#' @examples
#' df <- data.frame(
#'   assembly = LETTERS[1:5],
#'   continent = c("europe", "europe", "asia", "europe", "asia"),
#'   country = c("hungary", "hungary", "china", "serbia", "china"),
#'   mlst = c("ST1", "ST2", "ST1", "ST2", "ST3")
#' )
#' # which samples come from the same country
#' varid_matrix(df, var = "country")
#' # same but focus on European samples
#' varid_matrix(
#'   df, var = "country", focus_by = "continent", focus_on = "europe")
#' # which samples have the same mlst
#' varid_matrix(df, var = "mlst")
varid_matrix <- function(df,
                         var = "country",
                         focus_on = NULL,
                         focus_by = NULL) {
  if (length(var) != 1) {
    stop("Required argument 'var' must have length 1.")
  }
  if (var %in% names(df) == FALSE) {
    stop(paste0("Variable '", var, "' is missing from df."))
  }
  if (!is.null(focus_by)) {
    focus_by <- match.arg(focus_by, choices = names(df))
  }
  if (!is.null(focus_on)) {
    focus_on <- match.arg(focus_on, choices = unique(df[[focus_by]]))
  }
  varid <- matrix(0, nrow(df), nrow(df))
  for (i in unique(df[[var]])){
    index <- which(df[[var]] == i)
    varid[index, index] <- 1
  }
  if (!is.null(focus_on) & !is.null(focus_by)) {
    maskmat <- mask_matrix(df = df, focus_on = focus_on, focus_by = focus_by)
    varid <- varid * maskmat
  }
  diag(varid)<-NA
  rownames(varid) <- df$assembly
  colnames(varid) <- df$assembly
  return(varid)
}
