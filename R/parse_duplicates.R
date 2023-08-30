#' Parse names for identical sequences
#' 
#' This is a utility function for parsing the names of identical sequences from
#' a text file which was populated by seqkit.This function can be used for
#' removing identical tips from a phylogenetic tree before operations that are
#' sensitive to zero branch lengths and adding them back after the operation.
#' @param duplicates character; a character vector.
#' @return a list of duplicates. Each list element is a vector of ti
#' @examples
#' \dontrun{
#' # create an empty text file for storing duplicate sequence names
#' file.create(duplicates.txt)
#' # populate text file
#' seqkit rmdup -s seqs.fasta -D duplicates.txt -o seqs.nodup.fasta
#' # import duplicates
#' duplicates <- parse_duplicates("duplicates.txt")
#' }
parse_duplicates <- function(duplicates) {
  duplicates <- readLines(duplicates)
  if (length(duplicates) > 0) {
    duplicates <- gsub("^[0-9]+\\\t", "", duplicates) 
    duplicates <- strsplit(duplicates, ", ")
    names(duplicates) <- sapply(duplicates, function(x) x[1])
  } else {
    duplicates <- list()
  }
  return(duplicates)
}
