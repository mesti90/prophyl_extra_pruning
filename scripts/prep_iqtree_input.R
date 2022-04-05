rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]

input_dir <- paste0("./jobs/", jobname, "/output/gubbins/")
output_dir <- paste0("./jobs/", jobname, "/output/iqtree/")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

filenames <- c(
  "consensus.subs.node_labelled.final_tree.tre",
  "consensus.subs.filtered_polymorphic_sites.fasta"
)

for (i in filenames) {
  file.copy(from = paste0(input_dir, i), to = output_dir)
}
