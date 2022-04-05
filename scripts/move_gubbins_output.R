rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]

output_dir <- paste0("./jobs/", jobname, "/output/gubbins/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files <- list.files(full.names = TRUE)
files_sub <- files[grep("consensus", files)]

for (i in files_sub) {
  file.copy(from = i, to = output_dir)
  file.remove(i)
}
