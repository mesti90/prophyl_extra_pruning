rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]

output_dir <- paste0("./jobs/", jobname, "/output/gubbins")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

genomes <- list()
d1 <- list.dirs(output_dir, recursive = FALSE)
for (i in d1) {
  d2 <- list.dirs(i, recursive = FALSE)
  if (length(d2) > 0) {
    for (j in d2) {
      f <- seqinr::read.fasta(
        file = paste0(j, "/snps.consensus.subs.fa"),
        forceDNAtolower = FALSE)
      newname <- basename(j)
      #only keep chromosome - add later as decision point??
      genomes[[newname]] <- f[[1]]
    }
  }
}

seqinr::write.fasta(
  sequences = genomes,
  names = names(genomes),
  file.out = paste0(output_dir, "/consensus.subs.fasta")
)
