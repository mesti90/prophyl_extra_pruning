rm(list=ls())
library(ggtree) # required for tibble to convert tree to tibble

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
treemeta_file <- args[2]
ans_file <- args[3]

# load tree
tree <- ape::read.tree(tree_file)

# load metadata
meta <- read.csv(treemeta_file, sep = "\t")
meta <- dplyr::rename(meta, label = assembly)

# only keep metadata for assemblies that are on the tree
meta <- meta[which(meta$label %in% tree$tip.label),]

# load ancestral states
ans <- read.csv(ans_file, sep = "\t")
ans <- reshape2::dcast(ans, label ~ group)

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
  tree_tbl <- cbind(
    tree_tbl[,1:index],
    scp_df,
    tree_tbl[,(index+1):ncol(tree_tbl)]
  )
  
}


tree_tbl <- tibble::as_tibble(tree_tbl)
class(tree_tbl) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")




index_transmission_nodes <- vector()
for (i in 1:nrow(tree_tbl)) {
  country_child <- tree_tbl$country[i]
  country_parent <- tree_tbl$country[which(tree_tbl$node == tree_tbl$parent[i])]
  if (country_child == country_parent) next() else {
    index_transmission_nodes <- c(index_transmission_nodes, i)
  }
}

#tree_tbl$country_parent <- NA
#for (i in 1:nrow(tree_tbl)) {
#  index <- which(tree_tbl$node == tree_tbl$parent[i])
#  tree_tbl$country_parent[i] <- tree_tbl$country[index]
#}

tree_tbl$import <- FALSE
tree_tbl$import[index_transmission_nodes] <- TRUE


tree_tbl$collection_day_d <- lubridate::decimal_date(tree_tbl$collection_day)

mu <- ape::estimate.mu(tree, tree_tbl$collection_day_d)
dates <- ape::estimate.dates(tree, tree_tbl$collection_day_d, mu = mu)

tree_tbl$collection_day_pred <- as.Date(
  lubridate::date_decimal(dates),
  format = "%Y-%m-%d",
  origin = "1970-01-01"
)

tree_tbl$collection_year_pred <- lubridate::year(tree_tbl$collection_day_pred)

saveRDS(tree_tbl, file = "tree_tbl.rds")
