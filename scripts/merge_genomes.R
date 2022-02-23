rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]

output_dir <- paste0("./jobs/", jobname, "/output/snippy/")
genomenames <- list.dirs(output_dir, full.names= FALSE, recursive = FALSE)

genomes <- list()
for (i in genomenames) {
  f <- seqinr::read.fasta(
    file = paste0(output_dir, i, "/snps.consensus.subs.fa"),
    forceDNAtolower = FALSE)
  newname <- strsplit(i, "_")[[1]][1]
  
  #only keep chromosome - add later as decision point??
  genomes[[newname]] <- f[[1]]
}

seqinr::write.fasta(
  sequences = genomes,
  names = names(genomes),
  file.out = paste0(output_dir, "consensus.subs.fasta")
  )
