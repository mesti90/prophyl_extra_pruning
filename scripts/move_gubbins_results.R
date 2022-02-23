rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]

files <- list.files(full.names = TRUE)
files_sub <- files[grep("consensus", files)]

for (i in files_sub) {
  file.copy(from = i,
            to = paste0("./jobs/", jobname, "/output/gubbins/")) && file.remove(i)
}
