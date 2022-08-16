rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)

.libPaths(new = args[3])
snps = seqinr::read.fasta(args[1])
subsample_path = args[2]

subsample_index = strsplit(basename(subsample_path), "_|\\.")[[1]][2]
assemblies <- read.csv(subsample_path, sep = "\t")

index <- which(names(snps) %in% assemblies$assembly)

f <- snps[index]

seqinr::write.fasta(
  sequences = f,
  names = names(f),
  file.out = paste0("subsample_", subsample_index,".fasta")
)
