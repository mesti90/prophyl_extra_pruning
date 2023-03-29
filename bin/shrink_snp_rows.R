# Shrinking a phylogenetic tree may remove a number of tips from the tree. If
# this happens, the number of tips on the phylogenetic tree will be lower than
# the number of sequences in the fasta file. This script removes sequences from
# the fasta file to ensure all sequences match a tip of the phylogenetic tree.

rm(list = ls())

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  tree_path <- args[1]
  snps_path <- args[2]
} else {
  test_dir <- "~/Methods/prophyl-tests/test-shrink_snp_rows"
  tree_path <- paste0(test_dir, "/treeshrink.tre")
  snps_path <- paste0(test_dir, "/chromosomes.filtered_polymorphic_sites.fasta")
}

tree <- ape::read.tree(tree_path)
snps <- seqinr::read.fasta(snps_path, forceDNAtolower = FALSE)

index <- which(names(snps) %in% tree$tip.label == FALSE)

if (length(index) > 0) {
  snps <- snps[-index]
}

if (length(tree$tip.label) != length(snps)) {
  stop("Number of tips on tree is different from number of sequences in SNPs.")
}
if (any(tree$tip.label %in% names(snps) == FALSE)) {
  index <- which(tree$tip.label %in% names(snps) == FALSE)
  tips_collapsed <- paste(tree$tip.label[index], collapse = ", ")
  msg <- paste0("Some tips cannot be found within SNPs: ", tips_collapsed, ".")
  stop(msg)
}
if (any(names(snps) %in% tree$tip.label == FALSE)) {
  index <- which(names(snps) %in% tree$tip.label == FALSE)
  seqs_collapsed <- paste(names(snps)[index], collapse = ", ")
  msg <- paste0("Some sequences cannot be found on tree: ", seqs_collapsed, ".")
  stop()
}

if (!interactive()) {
  seqinr::write.fasta(
    sequences = snps,
    names = names(snps),
    file.out = "shrinked_snp_rows.fasta"
  )
}
