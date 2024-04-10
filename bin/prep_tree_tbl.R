library(ggtree) # required for tibble to convert tree to tibble
library(optparse)
rm(list = ls())

args <- commandArgs(trailingOnly = TRUE)

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "A tab delimited file with assembly metadata."
  ),
  make_option(
    c("-t", "--tree"),
    type = "character",
    help = "A tree in Newick format."
  ),
  make_option(
    c("-A", "--ancestral_states"),
    type = "character",
    help = "A tab delimited file with ancestral states."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/prophyl",
    assemblies = "assemblies.tsv",
    tree = "dated_tree.nwk",
    ancestral_states = "ancestral_states.tsv"
  )
}

library(devtools)
load_all(args$project_dir)

tree_file <- args$tree
treemeta_file <- args$assemblies
ans_file <- args$ancestral_states

# load tree
tree <- ape::read.tree(tree_file)

# load metadata
meta <- read.csv(treemeta_file, sep = "\t")
meta <- dplyr::rename(meta, label = assembly)

# only keep metadata for assemblies that are on the tree
meta <- meta[which(meta$label %in% tree$tip.label),]

# load ancestral states
ans <- read.csv(ans_file, sep = "\t")
ans <- tidyr::pivot_wider(ans, names_from = group, values_from = value)

# join tree and metadata
tree_tbl <- tibble::as_tibble(tree)
tree_tbl <- dplyr::left_join(tree_tbl, meta, by = "label")

# where ancestral states were predicted, replace columns with predicted
tree_tbl_names <- names(tree_tbl)
index <- which(names(tree_tbl) %in% names(ans))[-1]
tree_tbl <- tree_tbl[,-index]
tree_tbl <- dplyr::left_join(tree_tbl, ans, by = "label")
index <- unname(sapply(tree_tbl_names, function(x) which(names(tree_tbl) == x)))
tree_tbl <- tree_tbl[, index]

# combinatoric probability of state change
scp <- function(from, to){
  from_all <- strsplit(from, split = "\\|")[[1]]
  to_all <- strsplit(to, split = "\\|")[[1]]
  m <- matrix(0, nrow = length(from_all), ncol = length(to_all))
  rownames(m) <- from_all
  colnames(m) <- to_all
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      m[i,j] <- from_all[i] != to_all[j]
    }
  }
  p <- round(sum(m)/(length(from_all)*length(to_all)), 3)
  return(p)
}

for (i in names(ans)[-1]){
  parent <- NA
  for (j in 1:nrow(tree_tbl)) {
    index <- which(tree_tbl$node == tree_tbl$parent[j])
    parent[j] <- tree_tbl[[i]][index]
  }
  scp_result <- unname(mapply(function(x,y) scp(x,y), parent, tree_tbl[[i]]))
  scp_df <- data.frame(
    A = parent,
    B = scp_result
  )
  names(scp_df) <- c(paste0(i, "_from"), paste0(i, "_sc_prob"))
  index <- which(names(tree_tbl) == i)
  tf <- index < ncol(tree_tbl)
  if (tf) {
    tree_tbl <- cbind(
      tree_tbl[,1:index],
      scp_df,
      tree_tbl[,(index+1):ncol(tree_tbl)]
    )
  } else {
    tree_tbl <- cbind(tree_tbl, scp_df)
  }
}

tree_tbl$collection_day <- as.Date(tree_tbl$collection_day)

tree_tbl$collection_day <- ape::estimate.dates(
  tree, 
  node.dates = lubridate::decimal_date(tree_tbl$collection_day),
  mu = 1
)

tree_tbl$collection_day <- as.Date(
  lubridate::date_decimal(tree_tbl$collection_day), format = "%Y-%m-%d")

tree_tbl$collection_year <- as.numeric(
  substr(as.character(tree_tbl$collection_date), 1, 4))

tree_tbl <- tibble::as_tibble(tree_tbl)
class(tree_tbl) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")

write.table(
  tree_tbl,
  file = "tree_tbl.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

saveRDS(tree_tbl, file = "tree_tbl.rds")

write.table(
  tree_tbl,
  file = "tree_tbl.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)