#' Tidy duplicates for assembly subsets
#' 
#' When subsampling an assembly data set that may contain assemblies with
#' identical sequences, it is important to eliminate any duplicates and ensure
#' that the duplicates are represented by their reference assembly, i.e. the
#' assembly that is present in the sequence data set.
#' @param df data.frame; data set containing all samples
#' @param subset data.frame; a subset of the data set
#' @param id_var character; a variable in df and subset which contains the
#' duplicate names.
#' @param duplicates list; a list of duplicates. See \code{parse_duplicates()}
#' for more information.

#' @return a list with two elements, \code{"subset"} and \code{"duplist"}.
#' \code{"subset"} contains the tidied subset data frame where any duplicates
#' are replaced by their reference assembly and \code{"duplist"} contains the
#' assemblies that have been removed in the process.
tidy_duplicates <- function(df, subset, id_var, duplicates) {
  duplist <- list()
  if (length(duplicates) > 0) {
    names_duplist <- vector()
    k = 1
    # each element contains a vector of assemblies that are in the subset
    # they can be reference assemblies or not.
    for (i in seq_along(duplicates)) {
      index <- which(subset[[id_var]] %in% duplicates[[i]])
      if (length(index) > 0) {
        if (length(index) == 1 &&
            subset[[id_var]][index] == names(duplicates)[i]) {
          next()
        }
        duplist[[k]] <- subset[[id_var]][index]
        k = k + 1
      }
    }
    if (length(duplist) > 0) {
      # The name of each element is the name of the reference assembly.
      names(duplist) <- unname(sapply(duplist, function(x) {
        for (i in seq_along(duplicates)) {
          if (any(x %in% duplicates[[i]])) {
            return(names(duplicates)[i])
          }
        }
      }))
      # rearrange duplist elements to ensure if ref is included it is first
      for (i in seq_along(duplist)) {
        if (names(duplist)[i] %in% duplist[[i]]) {
          index <- which(duplist[[i]] == names(duplist[i]))
          duplist[[i]] <- c(duplist[[i]][index], duplist[[i]][-index])
        }
      }
      # Remove all assemblies that are listed as duplicates. This will remove
      # any references as well. At the same time add all references.
      index_remove <- unname(unlist(lapply(duplist, function(x) {
        which(subset$assembly %in% x) 
      })))
      index_add <- which(df[[id_var]] %in% names(duplist))
      subset <- dplyr::bind_rows(
        subset[-index_remove, ],
        df[index_add, ]
      )
    }
  }
  return(list("subset" = subset, "duplist" = duplist))
}
