#' Probability of state change
#' 
#' This functions calculated the combinatoric probability of state change
#' between two character strings.
#' @param from character; origin
#' @param to character; destination
#' @param from_which character
#' @param to_which character
#' @param split character; the string separating multiple states within
#' \code{from} and/or \code{to}.
#' @example 
#' scp("KL2|KL3|KL49", "KL2|KL3")
#' scp("KL2|KL3|KL49", "KL2|KL3", from_which = "KL2", to_which = "KL3")
scp <- function(from,
                to,
                from_which = NULL,
                to_which = NULL,
                split = "\\|"){
  foo <- function(arg1, arg2, arg3, arg4) {
    from_all <- strsplit(arg1, split = split)[[1]]
    to_all <- strsplit(arg2, split = split)[[1]]
    m <- matrix(0, nrow = length(from_all), ncol = length(to_all))
    rownames(m) <- from_all
    colnames(m) <- to_all
    for (i in 1:nrow(m)) {
      for (j in 1:ncol(m)) {
        m[i,j] <- from_all[i] != to_all[j]
      }
    }
    if (arg3 == "NULL") {
      index_from <- 1:nrow(m)
    } else {
      index_from <- which(rownames(m) %in% arg3)
    }
    if(arg4 == "NULL") {
      index_to <- 1:ncol(m)
    } else {
      index_to <- which(colnames(m) %in% arg4)
    }
    m_hit <- m[index_from, index_to]
    if("matrix" %in% class(m_hit) == FALSE) {
      m_hit <- as.matrix(m_hit, nrow = length(index_from), ncol = length(index_to))
      rownames(m_hit) <- rownames(m)[index_from]
      colnames(m_hit) <- colnames(m)[index_to]
    }
    p <- round(sum(m_hit)/(length(from_all)*length(to_all)), 3)
    return(p)
  }
  from_which <- ifelse(is.null(from_which), "NULL", from_which)
  to_which <- ifelse(is.null(to_which), "NULL", to_which)
  res <- mapply(foo, from, to, from_which, to_which)
  out <- unname(unlist(res))
  return(out)
}
