rm(list = ls())

library(optparse)

args_list <- list(
  make_option(
    c("-a", "--assembly_file"),
    type = "character",
    help = "Path to a tbl of assemblies"
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    assembly_file = "assemblies.tsv",
  )
}

# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

# read in assemblies
assemblies <- read.csv(args$assembly_file, sep = "\t")

# only keep assemblies where city of origin is known
# if (!"city" %in% names(assemblies)) {
#     stop("No 'city' column in assembly file")
# }
#
# index_drop <- which(is.na(assemblies$city))
# index_keep <- which(!is.na(assemblies$city))
#
# if (length(index_keep) == 0) {
#     stop("No assemblies with known city")
# } else {
#     message("Keeping ", length(index_keep), " assemblies with known city")
#     message("Dropping ", length(index_drop), " assemblies with unknown city")
# }
#
# assemblies <- assemblies[index_keep, ]

saveRDS(assemblies, file = "filtered_assemblies.rds")

# end logging
if (!interactive()) {
  sink(con)
}
