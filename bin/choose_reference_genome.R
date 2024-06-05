# This script chooses the reference genome for the analysis based on the
# longest contig length. If there are multiple genomes with the same longest
# contig length, the script chooses the genome with the earliest collection
# date. If there are multiple genomes with the same longest contig length and
# collection date, the script chooses the first genome in the list. This genome
# is then copied to the current directory and renamed to refgen_[genome_name].

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
load_all(args$project_dir)

df <- read_df(args$assemblies)

index_max_length <- which(df$longest_contig == max(df$longest_contig, na.rm = TRUE))
longest_assemblies <- df$assembly[index_max_length]

df_filtered <- df[which(df$assembly %in% longest_assemblies),]

if (nrow(df_filtered) > 1) {
  df_nona <- df_filtered[which(!is.na(df_filtered$collection_day)),]
  if (nrow(df_nona) == 0) {
    df_filtered <- df_filtered[1, ]
  }
  if (nrow(df_nona) == 1) {
    df_filtered <- df_nona
  }
  if (nrow(df_nona) > 1) {
    index <- which(df_nona$collection_day == min(df_nona$collection_day))
    df_filtered <- df_nona[index[1], ]
  }
}

file.copy(
  from = df_filtered$assembly_path,
  to = paste0("refgen_", basename(df_filtered$assembly_path))
)
