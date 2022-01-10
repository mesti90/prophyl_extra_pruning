rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

files <- dir(paste0("./input/", args[1], "/raw_reads/"))

index <- grep("R1", files)

pairs <- data.frame(
  "R1" = sort(files[index]),
  "R2" = sort(files[-index])
)

outdir <- paste0("./output/", args[1])
if (!dir.exists(outdir)) dir.create(outdir)

write.table(pairs,
            file = paste0(outdir, "/read_pairs.tsv"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
