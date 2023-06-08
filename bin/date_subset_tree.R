# This script takes
# - a subsample id
# - a phylogenetic tree with branch lengths in snp per site
# - an alignment used for calculating the phylogenetic tree
# - a table containing sample metadata
# - a value for the number of processors to use
# And calculates a dated phylogeny.

# If the collection date of a sample is known by day, it will be used as is. If
# the collection date of the sample is known by month or year, a matching range
# will be used as input and the middle value as the starting value. If the
# collection date of the sample is in another format or not known, the sample
# will be dropped from the analysis. This is different from the approach used
# by treedater authors in the article https://doi.org/10.1093/ve/vex025 where
# they seem to have interpolated unknown samples to fall between the boundaries
# of the rest of the genomes.

rm(list = ls())

library(dplyr)
library(magrittr)
library(treedater)

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  subset_id <- args[1]
  tree_path <- args[2]
  f_path <- args[3]
  assemblies_path <- args[4]
  ncpu <- as.numeric(args[5])
} else {
  subset_id <- "test"
  tree_path <- "rooted_subset_tree.rds"
  f_path <- "subsample_008_filtered.fasta"
  assemblies_path <- "assemblies.tsv"
  ncpu = 10
}

if (!interactive()) {
  # create log file and start logging
  if(!dir.exists(subset_id)) dir.create(subset_id)
  con <- file(paste0(subset_id, "/treedater_log.txt"))
  sink(con, split = TRUE)
}

# import subset tree.
# if rooting was based on rtt, this can be a list of trees
tree <- readRDS(tree_path)

# if the tree is NOT a list of trees, make it a list of length 1 so that the
# rest of the script can be used without modification
if ("tip.label" %in% names(tree)) {
  tree <- list(tree)
}

f <- seqinr::read.fasta(f_path)
assemblies <- read.csv(assemblies_path, sep = "\t", header = TRUE)

# drop tips which cannot be found in assembly table and give a warning
index <- which(tree[[1]]$tip.label %in% assemblies$assembly == FALSE)
if (length(index) > 0) {
  tips_to_drop <- tree[[1]]$tip.label[index]
  tips_to_drop_collapsed <- paste(tips_to_drop, collapse = ", ")
  for (i in 1:length(tree)) {
    tree[[i]] <- ape::drop.tip(tree[[i]], tips_to_drop)
  }
  msg <- paste0(
    "One or more tips could not be found in assembly table and were dropped: ",
    tips_to_drop_collapsed,
    "."
  )
  warning(msg)
}

# filter to assemblies that are included in the tree
index <- which(assemblies$assembly %in% tree[[1]]$tip.label == FALSE)
if (length(index) > 0) {
  assemblies <- assemblies[-index, ]
}

# filter to relevant columns
assemblies <- assemblies[, which(names(assemblies) %in% c(
  "assembly", "collection_date", "collection_day"
))]

date_lower <- function(dates) {
  foo <- function(x) {
    if (is.na(x)) return(NA)
    year = suppressWarnings(
      as.numeric(stringi::stri_sub(x, from = 1, to = 4))
    )
    if (is.na(year)) {
      return(NA)
    }
    date_elements <- strsplit(x, split = "-")[[1]]
    year <- date_elements[1]
    year_start <- as.Date(paste0(year, "-01-01"))
    year_end <- as.Date(paste0(year, "-12-31"))
    if (length(date_elements) == 3) {
      date <- as.Date(x, format = "%Y-%m-%d")
    }
    if (length(date_elements) == 2) {
      month <- date_elements[2]
      date <- as.Date(paste(year, month, "01", sep = "-"))
    }
    if (length(date_elements) == 1) {
      date <- as.Date(paste(year, "01", "01", sep = "-"))
    }
    date_out <- as.numeric(date-year_start)/as.numeric(year_end-year_start+1)
    # 1st january may be problematic
    date_out <- round(1000*date_out, 0)
    if(nchar(date_out) == 2){
      date_out <- paste0("0", date_out)
    }
    if(nchar(date_out) == 1) {
      date_out <- paste0("00", date_out)
    }
    date_out <- paste0(year, ".", date_out)
    
    return(date_out)
  }
  out <- unname(sapply(dates, function(x) try(foo(x), silent = TRUE)))
  out <- as.numeric(out)
  return(out)
}

date_upper <- function(dates) {
  foo <- function(x) {
    if(is.na(x)) return(NA)
    year = suppressWarnings(
      as.numeric(stringi::stri_sub(x, from = 1, to = 4))
    )
    if (is.na(year)) {
      return(NA)
    }
    date_elements <- strsplit(x, split = "-")[[1]]
    year <- date_elements[1]
    year_start <- as.Date(paste0(year, "-01-01"))
    year_end <- as.Date(paste0(year, "-12-31"))
    if (length(date_elements) == 3) {
      date <- as.Date(x, format = "%Y-%m-%d")
    }
    if (length(date_elements) == 2) {
      month <- date_elements[2]
      if (month == "02"){
        leap <- lubridate::leap_year(paste0(year, "-01-01"))
        if (leap) day <- "29"
        else day <- "28"
      }
      if (month != "02"){
        if (month %in% c("04", "06", "09", "11")) day <- "30"
        else day <- "31"
      }
      date <- as.Date(paste(year, month, day, sep = "-"))
    }
    if (length(date_elements) == 1) {
      date <- as.Date(paste(year, "12", "31", sep = "-"))
    }
    
    date_out <- as.numeric(date-year_start)/as.numeric(year_end-year_start+1)
    # 1st january may be problematic
    date_out <- round(1000*date_out, 0)
    if(nchar(date_out) == 2){
      date_out <- paste0("0", date_out)
    }
    if(nchar(date_out) == 1) {
      date_out <- paste0("00", date_out)
    }
    date_out <- paste0(year, ".", date_out)
    
    return(date_out)
  }
  out <- unname(sapply(dates, function(x) try(foo(x), silent = TRUE)))
  out <- as.numeric(out)
  return(out)
}

date_middle <- function(dates){
  foo <- function(x){
    if(is.na(x)) return(NA)
    year = suppressWarnings(
      as.numeric(stringi::stri_sub(x, from = 1, to = 4))
    )
    if (is.na(year)) {
      return(NA)
    }
    date_elements <- strsplit(x, split = "-")[[1]]
    year <- date_elements[1]
    if (length(date_elements) == 1){
      date <- paste(year, "06","15", sep = "-")
    }
    if (length(date_elements) == 2){
      month <- date_elements[2]
      date <- paste(year, month, "15", sep = "-")
    }
    if (length(date_elements) == 3){
      date <- x
    }
    return(date)
  }
  out <- as.Date(sapply(dates, foo))
  return(out)
}

# add new temporal variables
assemblies$date <- date_middle(assemblies$collection_date)
assemblies$lower <- date_lower(assemblies$collection_date)
assemblies$upper <- date_upper(assemblies$collection_date)

# define range for dates which are not known exactly
index_uncertain <- which(is.na(assemblies$collection_day))
uncertain_dates <- assemblies[index_uncertain, c("lower", "upper")]
rownames(uncertain_dates) <- assemblies$assembly[index_uncertain]

# drop tips where a range cannot be defined
index <- which(is.na(uncertain_dates$lower) | is.na(uncertain_dates$upper))
tips_to_drop <- rownames(uncertain_dates)[index]
if (length(index) > 0) {
  # drop from tree
  for (i in 1:length(tree)) {
    tree[[i]] <- ape::drop.tip(tree[[i]], tips_to_drop)
  }
  # drop from assemblies
  assemblies <- assemblies[-which(assemblies$assembly %in% tips_to_drop), ]
  # drop from uncertain dates
  uncertain_dates <- uncertain_dates[-which(row.names(uncertain_dates) %in% tips_to_drop), ]
  tips_to_drop_collapsed <- paste(tips_to_drop, collapse = ", ")
  msg <- paste0(
    "One or more tips were dropped because no data on sampling data was found: ",
    tips_to_drop_collapsed,
    "."
  )
  warning(msg)
}

# rename tips to include dates
for (i in 1:length(tree)) {
  tree[[i]]$tip.label <- sapply(tree[[i]]$tip.label, function(x) {
  index <- which(assemblies$assembly == x)
  paste(x, assemblies$date[index], sep = "|")
}, USE.NAMES = FALSE)
}

# rename rownames in uncertain dates to include dates
row.names(uncertain_dates) <- sapply(row.names(uncertain_dates), function(x) {
  index <- which(assemblies$assembly == x)
  paste(x, assemblies$date[index], sep = "|")
}, USE.NAMES = FALSE)

# extract dates from tip labels in appropriate format
sts <- sampleYearsFromLabels(tree[[1]]$tip.label, delimiter = "|")

# date tree

# if there are multiple rooted trees, date each one separately and then keep
# the one with the highest log likelihood
dtr <- list()
for (i in 1:length(tree)) {
  dtr[[i]] <- dater(
    tree[[i]],
    sts,
    s = length(f[[1]]),
    estimateSampleTimes = uncertain_dates,
    clock = 'strict', 
    ncpu =  ncpu
  )
}
loglik <- sapply(1:length(dtr), function(x) dtr[[x]]$loglik)
index <- which(loglik == max(loglik))[1]
dtr <- dtr[[index]]

# rescale non-dated branch lengths from per site to per genome (more intuitive)
alignment_length <- length(f[[1]])
dtr$intree$edge.length <- dtr$intree$edge.length*alignment_length

# export plots
try(dev.off(), silent = TRUE)
try(dev.off(), silent = TRUE)

pdf(file = paste0(subset_id, "/treedater_root_to_tip.pdf"))
rootToTipRegressionPlot(dtr)
dev.off()

try(dev.off(), silent = TRUE)
try(dev.off(), silent = TRUE)

png(file = paste0(subset_id, "/treedater_root_to_tip.png"))
rootToTipRegressionPlot(dtr)
dev.off()

# rename internal nodes
dtr %<>% makeNodeLabel(., method = "number", prefix = "Node_")

saveRDS(dtr, file = paste0(subset_id, "/dated_tree.rds"))

# do not remove dates from tip labels
# if dates are removed parboot will fail in a subsequent step
# dtr$tip.label <- unname(sapply(dtr$tip.label, function(x) {
#   strsplit(x, "\\|")[[1]][1]
# }))

ape::write.tree(dtr, file = paste0(subset_id, "/treedater_tree_with_time.nwk"))

if (!interactive()) {
  # end logging
  sink(con)
}
