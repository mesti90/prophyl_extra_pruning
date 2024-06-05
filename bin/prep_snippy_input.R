rm(list=ls())

if (!interactive()) {
  library(optparse)
  args_list <- list(
    make_option(
      "--project_dir",
      type = "character",
      help = "Path to the project directory",
      default = "prophyl-priv"
    ),
    make_option(
      "--assemblies",
      type = "character",
      help = "Path to a tbl of assemblies",
      default = "assemblies.tsv"
    )
  )
  args_parser  <- OptionParser(option_list = args_list)
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "prophyl-priv",
    assemblies = "assemblies.tsv"
  )
}

library(devtools)
library(dplyr)
load_all(args$project_dir)

df <- read_df(args$assemblies)

varnames <- c("assembly", "R1_path", "R2_path", "assembly_path")
index <- which(!varnames %in% colnames(df))

if (length(index) > 0) {
  varnames_missing_collapsed <- paste(varnames[index], collapse = ", ")
  msg <- paste0(
    "The following required variables are missing from the input table: ",
    varnames_missing_collapsed,
    "."
  )
  stop(msg)
}

paired_reads <- data.frame()
single_reads <- data.frame()
contigs <- data.frame()
error <- data.frame()

for (i in 1:nrow(df)) {
  if (!is.na(df$R1_path[i])) {
    if (!is.na(df$R2_path)) {
      # append paired_reads
      R1_exists <- file.exists(df$R1_path[i])
      R2_exists <- file.exists(df$R2_path[i])
      if (all(R1_exists, R2_exists)) {
        new_row <- data.frame(
          assembly = df$assembly[i],
          R1 = df$R1_path[i],
          R2 = df$R2_path[i]
        )
        paired_reads <- dplyr::bind_rows(paired_reads, new_row)
      } else {
        new_row <- data.frame(
          assembly = df$assembly[i],
          reason = "R1 or R2 file not found."
        )
        error <- dplyr::bind_rows(error, new_row)
      }
    } else {
      # append single_reads
      R1_exists <- file.exists(df$R1_path[i])
      if (R1_exists) {
        new_row <- data.frame(
          assembly = df$assembly[i],
          reads = df$R1_path[i]
        )
        single_reads <- dplyr::bind_rows(single_reads, new_row)
      } else {
        new_row <- data.frame(
          assembly = df$assembly[i],
          reason = "R1 file not found."
        )
        error <- dplyr::bind_rows(error, new_row)
      }
    }
  } else if (!is.na(df$R2_path)) {
    new_row <- data.frame(
      assembly = df$assembly[i],
      reason = "Syntax error. For single end reads use R1."
    )
    error <- dplyr::bind_rows(error, new_row)
  } else if (!is.na(df$assembly_path[i])) {
    # append contigs
    assembly_exists <- file.exists(df$assembly_path[i])
    if (assembly_exists) {
      new_row <- data.frame(
        assembly = df$assembly[i],
        path = df$assembly_path[i]
      )
      contigs <- dplyr::bind_rows(contigs, new_row)
    } else {
      new_row <- data.frame(
        assembly = df$assembly[i],
        reason = "Assembly file not found."
      )
      error <- dplyr::bind_rows(error, new_row)
    }
  } else {
    new_row <- data.frame(
      assembly = df$assembly[i],
      reason = "No read or assembly files found."
    )
    error <- dplyr::bind_rows(error, new_row)
  }
}

write_tsv(paired_reads, "paired_reads.tsv")
write_tsv(single_reads, "single_reads.tsv")
write_tsv(contigs, "contigs.tsv")
write_tsv(error, "error.tsv")

if (nrow(error) > 0) {
  stop("Some files were not found. See error.tsv for details.")
}
