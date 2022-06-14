rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
assembly_name <- args[1]
assembly_dir <- args[2]

genomes <- list()
f <- seqinr::read.fasta(
  file = paste0(assembly_dir, "/snps.consensus.subs.fa"),
  forceDNAtolower = FALSE)
newname <- basename(assembly_name)
genomes[[newname]] <- f[[1]]

seqinr::write.fasta(
  sequences = genomes,
  names = names(genomes),
  file.out = paste0(assembly_name, "_wgs.fasta")
)
