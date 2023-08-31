validate_focus <- function(df, focus_by, focus_on) {
  # validate focus_by and focus_on together
  if (sum(is.null(focus_on), is.null(focus_by)) == 1) {
    stop("Either provide both 'focus_by' and 'focus_on' or none.") 
  }
  # validate focus_by
  if (!is.null(focus_by)) {
    focus_by <- match.arg(focus_by, choices = names(df))
  }
  # validate focus_on
  if (!is.null(focus_on)) {
    focus_on <- match.arg(
      focus_on, choices = sort(unique(df[[focus_by]])), several.ok = TRUE)
  }
  return(list("focus_by" = focus_by, "focus_on" = focus_on))
}

validate_id_var <- function(df, id_var) {
  id_var <- match.arg(id_var, choices = names(df))
  if (length(df[[id_var]]) > length(unique(df[[id_var]]))) {
    msg <- paste0("Variable ", id_var," contains duplicates.")
    stop(msg)
  }
  if (any(is.na(df[[id_var]]))) {
    msg <- paste0("Variable ", id_var," contains missing values.")
    stop(msg)
  }
  return(id_var)
}