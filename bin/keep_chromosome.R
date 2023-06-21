rm(list=ls())
library(seqinr)

args <- commandArgs(trailingOnly = TRUE)

f <- seqinr::read.fasta(
  file = paste0(args[1], "/snps.consensus.subs.fa"),
  forceDNAtolower = FALSE
)

f_lengths <- sapply(f, length)

seqinr::write.fasta(
  sequences = f[[which.max(f_lengths)]],
  names = basename(args[1]),
  file.out = paste0(basename(args[1]), ".fasta")
)
