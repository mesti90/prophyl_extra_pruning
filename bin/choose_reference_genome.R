# This script chooses the reference genome for the analysis based on the
# longest contig length. If there are multiple genomes with the same longest
# contig length, the script chooses the genome with the earliest collection
# date. If there are multiple genomes with the same longest contig length and
# collection date, the script chooses the first genome in the list. This genome
# is then copied to the current directory and renamed to refgen_[genome_name].
# 
# Update: The earlier protocol above sometimes led to erroneous result, e.g.
# with klebsi ST383 where the selected reference genome was unrealistically 
# large for some reason and then after mapping almost all chromosomes turned 
# out to be identical which is clearly wrong. Instead I have added that if
# possible, we only look at GCF genomes, and we always choose a reference that
# is not extremely long.

rm(list=ls())

library(optparse)
args_list <- list(
  make_option(
    "--project_dir",
    type = "character",
    help = "Path to the project directory",
    default = "prophyl"
  ),
  make_option(
    "--assemblies",
    type = "character",
    help = "Path to a tbl of assemblies",
    default = "data/tree_input/tree_input_ST383.tsv"
  )
)

args_parser  <- OptionParser(option_list = args_list)
args  <- parse_args(args_parser)

library(devtools)
library(R.utils)
load_all(args$project_dir)

df <- read_df(args$assemblies)

# Filter to GCF genomes, if there are any
index_gcf <- grep("^GCF", df$assembly)
if (length(index_gcf) > 0) df <- df [index_gcf,]

# Filter to genomes that are not too long
largest_accepted_genome <- get_non_outlier_range(
  df$genome_size, qlow = 0, qhigh = 0.95, m = 1, keep = "lower"
)[2]
df <- df[which(df$genome_size <= largest_accepted_genome),]

index_max_length <- which(df$longest_contig == max(df$longest_contig, na.rm = TRUE))
longest_assemblies <- df$assembly[index_max_length]

df_filtered <- df[which(df$assembly %in% longest_assemblies),]

# collection_day
df_filtered$collection_day <- unname(date_middle(df_filtered$collection_date))

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

infile <- df_filtered$assembly_path
local_copy <- paste0("refgen_", basename(infile))

file.copy(from = infile, to = local_copy)

# if file is compressed, uncompress it

filetype = summary(file(local_copy))$class

if (filetype == "gzfile") {
  R.utils::gunzip(
    filename = local_copy,
    skip = TRUE,
    remove = TRUE
  )
}
