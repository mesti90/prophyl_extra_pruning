#' Create symbolic links
#' 
#' This function creates a directory structure within the working directory for
#' storing genome sequences in a standardized form, and then creates symbolic
#' links which point to genome sequences elsewhere in the file system.
#' @param from charcater; path to a directory which contains genome sequences.
#' @param type character; the type of genomes sequences stored in the
#' \code{from} directory. Can be either \code{"raw_reads"},
#' or \code{"assembled_genomes"}.
#' @param fsep character; the path separator to use.
#' @details The analysis pipeline requires that genome sequences are organised
#' in a standardized directory structure inside the working directory. This
#' structure is designed to provide maximum flexibility. While it is possible to
#' organize the genome sequences themselves, it is more efficient to create
#' symbolic links instead which point to these files. This function sets up the
#' directory structure inside the working directory and fills it up with
#' symbolic links.
#' @examples 
#' \dontrun{
#' # original directory contains raw_reads
#' create_link(
#'   from = "path_before_project_name/project_name/genome_dir",
#'   type = "raw_reads")
#' # original directory contains assembled genomes
#' create_link(
#'   from = "path_before_project_name/project_name/genome_dir",
#'   type = "assembled_genomes") 
#' }
create_link <- function(from, type, fsep = "/") {
  foo <- function(x, y) {
    y <- match.arg(y, c("raw_reads", "assembled_genomes"))
    path <- strsplit(x, split = fsep)[[1]]
    to <- paste0(getwd(), "/genomes/", path[length(path)-1], "/", y, "/")
    if (!dir.exists(to)) {
      dir.create(to, recursive = TRUE)
      if (grepl("/$", x) == FALSE) x <- paste0(x, "/")
      system(paste0("ln -s ", x, "*", " --target ", to))
    }
  }
  out <- mapply(foo, from, type)
}
