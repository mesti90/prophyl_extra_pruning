rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

.libPaths(new = args[2])
library(devtools)
library(lubridate)
load_all()

assemblies <- read.csv(args[1], sep = "\t")
assemblies <- assemblies[order(assemblies$assembly),]

# temporal distance - time difference between collections dates

dates <- unname(lubridate::decimal_date(date_middle(assemblies$collection_date)))
colldist = round(abs(outer(dates, dates, "-")),2)
colnames(colldist) <- assemblies$assembly
row.names(colldist) <- assemblies$assembly
saveRDS(colldist, file = "colldist.rds")