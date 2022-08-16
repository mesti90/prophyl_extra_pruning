rm(list=ls())
library(dplyr)
library(magrittr)
library(treedater)

args <- commandArgs(trailingOnly = TRUE)

subset_id <- args[1]
tree <- ape::read.tree(file = args[2])
f <- seqinr::read.fasta(file = args[3])
assemblies <- read.csv(args[4], sep = "\t", header = TRUE)
ncpu <- as.numeric(args[5])

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

if(!dir.exists(subset_id)) dir.create(subset_id)

con <- file(paste0(subset_id, "/treedater_log.txt"))
sink(con, split = TRUE)

dtr <- dater(tree,
             sts,
             s = length(f[[1]]),
             clock = 'strict', 
             ncpu =  ncpu)

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

dtr %<>% makeNodeLabel(., method = "number", prefix = "Node_")

saveRDS(dtr, file = paste0(subset_id, "/dated_tree.rds"))

dtr$tip.label <- unname(sapply(dtr$tip.label, function(x) {
  strsplit(x, "\\|")[[1]][1]
}))

ape::write.tree(dtr, file = paste0(subset_id, "/treedater_tree_with_time.nwk"))

sink(con)
