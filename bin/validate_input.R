rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
df <- read.csv(args[1], sep = "\t")

# if pastml would break because of too many character states, stop.
# look at how this error appears in pastml, copy.

if ("assembly" %in% names(df) == FALSE) {
  stop("Required column 'assembly' is missing.")
}

if (length(unique(df$assembly)) != length(df$assembly)) {
  stop("Assembly names must be unique. Check.")
}

if (sum(grepl("root", df$assembly, ignore.case = TRUE)) > 0) {
  stop("Assembly name cannot be 'root'. Specify another name.")
}

if (sum(grepl("^node_", df$assembly, ignore.case = TRUE)) > 0) {
  stop("Assembly name cannot start with 'Node_'. Specify another name.")
}

if (sum(!is.na(suppressWarnings(as.numeric(df$assembly)))) > 0) {
  stop("Assembly name cannot be a number. Specify another name.")
}

if ("collection_day" %in% names(df) == FALSE) {
  stop("Required column 'collection_day' is missing.")
}

#if (class(df$collection_day) != "Date") {
#  stop("Required column 'collection_day' must be a 'Date'.")
#}

# if (sum(is.na(df$collection_day)) > 0) {
#   stop("All assemblies must contain collection days. Check.")
# }

if ("country" %in% names(df) == FALSE) {
  stop("Required column 'country' is missing.")
}

if ("country_iso2c" %in% names(df) == FALSE) {
  stop("Required column 'country_iso2c' is missing.")
}

if ("continent" %in% names(df) == FALSE) {
  stop("Required column 'continent' is missing.")
}

# Required variables for ancestral state prediction
# TODO expand this structure to all of the script
required_variables <- c(
  "country", "region23", "continent", "k_serotype", "k_confidence"
)
index <- which(required_variables %in% names(df) == FALSE)
if (length(index) > 0) {
  missing = paste(required_variables[index], collapse = ", ")
  stop("One or mode required variables are missing: ", missing, ".")
}



write.table(
  df, 
  file = "assemblies.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
