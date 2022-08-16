rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

assemblies <- read.csv(args[1], sep = "\t")
nsubs <- as.numeric(args[2])
ss <- round(nrow(assemblies)*as.numeric(args[3]), 0)

set.seed(0)
digits <- ceiling(log10(nsubs))
for (i in 1:nsubs) {
  zeroes <- digits - floor(log10(i))
  write.table(
    assemblies[sample(1:nrow(assemblies), ss, replace = FALSE), ],
    file = paste0("subsample_", rep(0, times = zeroes), i, ".tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
