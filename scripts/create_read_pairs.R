rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]

files <- dir(paste0("./jobs/", jobname, "/raw_reads/"))

index <- grep("R1", files)

pairs <- data.frame(
  "R1" = sort(files[index]),
  "R2" = sort(files[-index])
)

output_dir <- paste0("./jobs/", jobname, "/output/")
if (!dir.exists(output_dir)) dir.create(output_dir)

write.table(pairs,
            file = paste0(output_dir, "read_pairs.tsv"),
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")
