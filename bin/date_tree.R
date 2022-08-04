rm(list=ls())
library(dplyr)
library(magrittr)
library(treedater)

args <- commandArgs(trailingOnly = TRUE)

tree <- ape::read.tree(file = args[1])
f <- seqinr::read.fasta(file = args[2])
assemblies <- read.csv(args[3], sep = "\t", header = TRUE)
ncpu <- as.numeric(args[4])

dates <- assemblies[, which(names(assemblies) %in% c(
  "assembly", "collection_day"
))]
dates <- dplyr::rename(dates, name = assembly)
dates <- dplyr::rename(dates, date = collection_day)

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

con <- file("treedater_log.txt")
sink(con, split = TRUE)

dtr <- dater(tree,
             sts,
             s = length(f[[1]]),
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
