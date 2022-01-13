rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

files <- list.files(full.names = TRUE)
files_sub <- files[grep("consensus", files)]

for (i in files_sub) {
  file.copy(from = i,
            to = paste0("./output/", args[1], "/gubbins/")) && file.remove(i)
}
