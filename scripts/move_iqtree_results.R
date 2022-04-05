rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]

files <- list.files(full.names = TRUE)

include <- c(
  "\\.contree$",
  "\\.iqtree$",
  "\\.log$",
  "\\.treefile$",
  "\\.ckp.gz$",
  "\\.model.gz",
  "\\.splits.nex",
  "\\.uniqueseq.phy"
)

files_sub <- unname(sapply(include, function(x) {
  files[grep(x, files)]
}))

for (i in files_sub) {
  file.copy(from = i,
            to = paste0("./jobs/", jobname, "/output/iqtree/")) && file.remove(i)
}
