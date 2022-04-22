rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

jobname <- args[1] # name of the directory inside the input directory
p_cores <- as.numeric(args[2]) # proportion of cores to use for the analysis

library(dplyr)
library(treedater)

# create output directory
tdir <- paste0("./jobs/",jobname,"/output/treedater/")
if (!dir.exists(tdir)) dir.create(tdir, recursive = TRUE)

# read tree from gubbins directory
gdir <- paste0("./jobs/",jobname,"/output/gubbins/")
tree <- ape::read.tree(
  file = paste0(gdir, "consensus.subs.node_labelled.final_tree.tre"))

# unroot tree, if rooted
if (ape::is.rooted(tree)) tree <- ape::unroot(tree)

# read data set with dates, requires full dates

dates <- read.table(paste0("./jobs/", jobname,"/dates.tsv"),
                    sep = '\t',
                    header = TRUE)
dates$date <- as.Date(dates$date)

# rename tips to include dates

tree$tip.label <- sapply(tree$tip.label, function(x) {
  index <- which(dates$name == x)
  paste(x, dates$date[index], sep = "|")
}, USE.NAMES = FALSE)

# extract dates from tip labels in appropriate format

sts <- sampleYearsFromLabels(tree$tip.label, delimiter = "|")

f <- seqinr::read.fasta(
  file = paste0(gdir, "consensus.subs.filtered_polymorphic_sites.fasta"))

con <- file(paste0(tdir, "treedater_log.txt"))
sink(con, split = TRUE)

dtr <- dater(tree,
             sts,
             s = length(f[[1]]),
             clock = 'strict', 
             ncpu =  max(round(p_cores*parallel::detectCores(),0),1))

try(dev.off(), silent = TRUE)
try(dev.off(), silent = TRUE)

pdf(file = paste0(tdir, "treedater_root_to_tip.pdf"))
rootToTipRegressionPlot(dtr)
dev.off()

try(dev.off(), silent = TRUE)
try(dev.off(), silent = TRUE)

png(file = paste0(tdir, "treedater_root_to_tip.png"))
rootToTipRegressionPlot(dtr)
dev.off()

dtr$tip.label <- unname(sapply(dtr$tip.label, function(x) {
  strsplit(x, "\\|")[[1]][1]
}))

ape::write.tree(dtr,
                file = paste0(tdir, "treedater_tree_with_time.nwk"))

sink(con)
