#' Probability of state change
#' 
#' This function calculates the combinatoric probability of state change
#' between two character strings.
#' @param from character; origin
#' @param to character; destination
#' @param from_which character
#' @param to_which character
#' @param split character; the string separating multiple states within
#' \code{from} and/or \code{to}.
#' @examples 
#' scp("KL2|KL3|KL49", "KL2|KL3")
#' scp("KL2|KL3|KL49", "KL2|KL3", from_which = "KL2", to_which = "KL3")
#' scp("KL2", "KL2|KL3")
#' @export
scp <- function(from,
                to,
                from_which = NULL,
                to_which = NULL,
                split = "\\|"){
  from_all <- strsplit(from, split = split)[[1]]
  to_all <- strsplit(to, split = split)[[1]]
  m <- matrix(0, nrow = length(from_all), ncol = length(to_all))
  rownames(m) <- from_all
  colnames(m) <- to_all
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      m[i,j] <- from_all[i] != to_all[j]
    }
  }
  if (is.null(from_which)) {
    index_from <- 1:nrow(m)
  } else {
    index_from <- which(rownames(m) %in% from_which)
  }
  if(is.null(to_which)) {
    index_to <- 1:ncol(m)
  } else {
    index_to <- which(colnames(m) %in% to_which)
  }
  m_hit <- m[index_from, index_to]
  if("matrix" %in% class(m_hit) == FALSE) {
    m_hit <- as.matrix(m_hit, nrow = length(index_from), ncol = length(index_to))
    if (length(index_to) > 1) {
      m_hit <- t(m_hit)
    }
    rownames(m_hit) <- rownames(m)[index_from]
    colnames(m_hit) <- colnames(m)[index_to]
  }
  p <- round(sum(m_hit)/(length(from_all)*length(to_all)), 3)
  out <- list(m = m, prob = p)
  return(out)
}
