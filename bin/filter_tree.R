rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)

tree_path <- args[1]
assembly_path <- args[2]

tree <- ape(read.tree(tree_path))
assemblies <- read.csv(assembly_path, sep = "\t")

# only keep assemblies where k_serotype is not NA
index <- which(is.na(assemblies$k_serotype))
if (length(index) > 0) {
    assemblies <- assemblies[-index, ]
}

# only keep assemblies where k_confidence is at least "Good"
index <- which(assemblies$k_confidence %in% c("None", "Low"))
if (length(index) > 0) {
    assemblies <- assemblies[-index, ]
}

tree <- ape::keep.tip(tree, assemblies$assembly)

ape::write.tree(tree, file = "filtered_tree.nwk")
saveRDS(assembly, file = "assembly_filtered.rds")
