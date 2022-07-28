rm(list=ls())
library(seqinr)

args <- commandArgs(trailingOnly = TRUE)

f <- seqinr::read.fasta(
  file = paste0(args[1], "/snps.consensus.subs.fa"),
  forceDNAtolower = FALSE
)

seqinr::write.fasta(
  sequences = f[[1]],
  names = basename(args[1]),
  file.out = paste0(basename(args[1]), ".fasta")
)
