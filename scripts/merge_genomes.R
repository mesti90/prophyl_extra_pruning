rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

outdir <- paste0("./output/", args[1], "/snippy")
genomenames <- list.dirs(outdir, full.names= FALSE, recursive = FALSE)

genomes <- list()
for (i in genomenames) {
  f <- seqinr::read.fasta(
    file = paste0(outdir, "/", i, "/snps.consensus.subs.fa"),
    forceDNAtolower = FALSE)
  newname <- strsplit(i, "_")[[1]][1]
  
  #only keep chromosome - add later as decision point??
  genomes[[newname]] <- f[[1]]
}

seqinr::write.fasta(
  sequences = genomes,
  names = names(genomes),
  file.out = paste0("./output/", args[1], "/snippy/consensus.subs.fasta")
  )
