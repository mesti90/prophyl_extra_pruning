#' Tidy duplicates for assembly subsets
#' 
#' When subsampling an assembly data set that may contain assemblies with
#' identical sequences, it is important to eliminate any duplicates and ensure
#' that the duplicates are represented by their reference assembly, i.e. the
#' assembly that is present in the sequence data set.
#' @param assemblies data.frame; assembly data set containing all samples
#' @param subset data.frame; a subset of the assembly data set
#' @param duplicates list; a list of duplicates. See \code{parse_duplicates()}
#' for more information.
#' @return a list with two elements, \code{"subset"} and \code{"duplist"}.
#' \code{"subset"} contains the tidied subset data frame where any duplicates
#' are replaced by their reference assembly and \code{"duplist"} contains the
#' assemblies that have been removed in the process.
tidy_duplicates <- function(assemblies, subset, duplicates) {
  duplist <- list()
  if (length(duplicates) > 0) {
    names_duplist <- vector()
    k = 1
    for (i in seq_along(duplicates)) {
      index <- which(subset$assembly %in% duplicates[[i]])
      if (length(index) > 0) {
        if (length(index) == 1 && 
            subset$assembly[index] == names(duplicates)[i]) {
          next()
        }
        duplist[[k]] <- subset$assembly[index]
        subset <- subset[-index, ]
        k = k + 1
      }
    }
    if (length(duplist) > 0) {
      names(duplist) <- unname(sapply(duplist, function(x) {
        for (i in seq_along(duplicates)) {
          if (any(x %in% duplicates[[i]])) {
            return(names(duplicates)[i])
          }
        }
      }))
      subset <- dplyr::bind_rows(
        subset,
        assemblies[which(assemblies$assembly %in% names(duplist)), ]
      )
      for (i in seq_along(duplist)) {
        if (names(duplist)[i] %in% duplist[[i]]) {
          index <- which(duplist[[i]] == names(duplist[i]))
          duplist[[i]] <- c(duplist[[i]][index], duplist[[i]][-index])
        }
      }
    }
  }
  return(list("subset" = subset, "duplist" = duplist))
}
