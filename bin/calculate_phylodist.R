rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

simtrees <- readRDS(args[1])
colldist <- readRDS(args[2])

nsim <- length(simtrees$trees)
ntips <- length(simtrees$trees[[1]]$tip.label)
phylodist_array <- array(NaN, c(ntips, ntips, nsim))

for (i in 1:nsim) {
  phylodist <- ape::cophenetic.phylo(simtrees$trees[[i]])
  index <- order(colnames(phylodist))
  phylodist <- phylodist[index, index]
  diag(phylodist) <- NA
  
  if (mean(colnames(phylodist) == colnames(colldist)) < 1) {
    stop("Distance matrix names do not match. Check.")
  }
  mrca = (phylodist - colldist)/2
  mrca[which(mrca < 0)] = 0
  phylodist_array[, , i] <- mrca
}

saveRDS(phylodist_array, file = "phylodist_array.rds")
