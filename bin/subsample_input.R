rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

assemblies <- read.csv(args[1], sep = "\t")
subsample_count <- as.numeric(args[2])
subsample_tipcount <- as.numeric(args[3])

type <- "balanced"
balance_by <- "region23"

if (balance_by %in% names(assemblies) == FALSE) {
  stop("Subsampling failed, variable '", balance_by, "' not found.")
}

if (type == "random") {
  set.seed(0)
  digits <- ceiling(log10(subsample_count+1))
  for (i in 1:subsample_count) {
    zeroes <- digits - floor(log10(i))
    filename = paste0(c("subsample_", rep(0, times = zeroes), i, ".tsv"), collapse = "")
    write.table(
      assemblies[sample(1:nrow(assemblies), subsample_tipcount, replace = FALSE), ],
      file = filename,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
  }
}

if (type == "balanced") {
  set.seed(0)
  years_min <- min(assemblies$collection_year, na.rm = TRUE)
  years_max <- max(assemblies$collection_year, na.rm = TRUE)
  years <- years_min:years_max
  
  size = 0
  k = 0
  
  while (size < subsample_tipcount & size < nrow(assemblies)) {
    k = k + 1
    index <- vector()
    for (i in years) {
      bins <- unique(assemblies[[balance_by]][which(assemblies$collection_year == i)])
      for (j in bins) {
        idx <- which(assemblies$collection_year == i & assemblies[[balance_by]] == j)
        if (length(idx) >= k) {
          index <- c(index, sample(idx, k))
        } else {
          index <- c(index, idx)
        }
      }
    }
    size <- length(index)
  }
  if (k == nrow(assemblies)) {
    stop("Subsampling did not converge. Check.")
  } else {
    digits <- ceiling(log10(subsample_count+1))
    for (h in 1:subsample_count) {
      index <- vector()
      for (i in years) {
        bins <- unique(assemblies[[balance_by]][which(assemblies$collection_year == i)])
        for (j in bins) {
          idx <- which(assemblies$collection_year == i & assemblies[[balance_by]] == j)
          if (length(idx) >= k) {
            index <- c(index, sample(idx, k))
          } else {
            index <- c(index, idx)
          }
        }
      }
      zeroes <- digits - floor(log10(h))
      filename = paste0(c("subsample_", rep(0, times = zeroes), h, ".tsv"), collapse = "")
      write.table(
        assemblies[index, ],
        file = filename,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  }
}
