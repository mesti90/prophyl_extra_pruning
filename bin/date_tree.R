rm(list=ls())
library(dplyr)
library(magrittr)
library(treedater)

args <- commandArgs(trailingOnly = TRUE)

tree <- ape::read.tree(file = args[1])
f <- seqinr::read.fasta(file = args[2])
assemblies <- read.csv(args[3], sep = "\t", header = TRUE)

ncpu <- as.numeric(args[4])

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

dates <- assemblies[, which(names(assemblies) %in% c(
  "assembly", "collection_date", "collection_day"
))]
dates <- dplyr::rename(dates, name = assembly)
dates$date <- date_middle(dates$collection_date)
dates$lower <- date_lower(dates$collection_date)
dates$upper <- date_upper(dates$collection_date)

# unroot tree, if rooted
if (ape::is.rooted(tree)) tree <- ape::unroot(tree)

# read data set with dates, requires full dates
dates$date <- as.Date(dates$date)

# rename tips to include dates

tree$tip.label <- sapply(tree$tip.label, function(x) {
  index <- which(dates$name == x)
  paste(x, dates$date[index], sep = "|")
}, USE.NAMES = FALSE)

# extract dates from tip labels in appropriate format

sts <- sampleYearsFromLabels(tree$tip.label, delimiter = "|")

# define range for dates which are not known exactly

index_uncertain <- which(is.na(dates$collection_day))
uncertain_dates <- dates[index_uncertain, c("lower", "upper")]
rownames(uncertain_dates) <- dates$name[index_uncertain]

# drop tips where a range cannot be defined

index <- which(is.na(uncertain_dates$lower) | is.na(uncertain_dates$upper))
tree <- ape::drop.tip(tree, rownames(uncertain_dates)[index])

# date tree

con <- file("treedater_log.txt")
sink(con, split = TRUE)

dtr <- dater(tree,
             sts,
             s = length(f[[1]]),
             estimateSampleTimes = uncertain_dates,
             clock = 'strict', 
             ncpu =  ncpu)

try(dev.off(), silent = TRUE)
try(dev.off(), silent = TRUE)

pdf(file = "treedater_root_to_tip.pdf")
rootToTipRegressionPlot(dtr)
dev.off()

try(dev.off(), silent = TRUE)
try(dev.off(), silent = TRUE)

png(file = "treedater_root_to_tip.png")
rootToTipRegressionPlot(dtr)
dev.off()

dtr %<>% makeNodeLabel(., method = "number", prefix = "Node_")

saveRDS(dtr, file = "dated_tree.rds")

dtr$tip.label <- unname(sapply(dtr$tip.label, function(x) {
  strsplit(x, "\\|")[[1]][1]
}))

ape::write.tree(dtr, file = "treedater_tree_with_time.nwk")

sink(con)
