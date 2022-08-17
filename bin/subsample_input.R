rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

assemblies <- read.csv(args[1], sep = "\t")
nsubs <- as.numeric(args[2])
ss <- as.numeric(args[3])

set.seed(0)
digits <- ceiling(log10(nsubs+1))
for (i in 1:nsubs) {
  zeroes <- digits - floor(log10(i))
  filename = paste0(c("subsample_", rep(0, times = zeroes), i, ".tsv"), collapse = "")
  write.table(
    assemblies[sample(1:nrow(assemblies), ss, replace = FALSE), ],
    file = filename,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
