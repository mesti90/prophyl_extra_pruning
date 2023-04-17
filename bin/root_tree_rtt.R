# This script takes a phylogenetic tree and attempts to root it using Minimum
# Ancestor Deviation (MAD).

rm(list = ls())

# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

# define input
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  project_dir <- args[1]
  tree_path <- args[2]
  assemblies_path <- args[3]
  threads <- args[4]
} else {
  test_dir <- "~/Methods/prophyl-tests/test-root_tree_rtt"
  project_dir <- "~/Methods/prophyl"
  tree_path <- paste0(test_dir, "/sample_tree.tre")
  assemblies_path <- paste0(test_dir, "/assemblies.tsv")
  threads <- 2
}

library(devtools)
library(dplyr)
library(ggplot2)
load_all(project_dir)

# read tree
tree <- ape::read.tree(tree_path)

# read assemblies
assemblies <- read.csv(assemblies_path, sep = "\t")

# ensure collection_day is a 'Date'
assemblies$collection_day <- date_middle(assemblies$collection_date)

# if collection_day is NA, remove both from assembly tbl and tree
remove <- assemblies$assembly[which(is.na(assemblies$collection_day))]

# remove from assembly tbl
assemblies <- assemblies[-which(assemblies$assembly %in% remove),]

# remove from tree
tree <- ape::drop.tip(tree, remove)

# collect tip dates in the same order as tree$tip.label
tip_dates <- unname(sapply(tree$tip.label, function(x) {
  index <- which(assemblies$assembly == x)
  assemblies$collection_day[index]
}))

# if tree is rooted, unroot
if (ape::is.rooted(tree)) {
  tree <- ape::unroot(tree)
}

# define custom objective function, perform robust linear regression
# TODO find another metric that works with residuals, is standardised, but
# penalises negative deviations more than positive deviations
objective_rlm <- function(x,y) -summary(MASS::rlm(y ~ x))$sigma^2

objective <- list(
  "correlation" = NULL,
  "rsquared" = NULL,
  "rms" = NULL,
  "rlm" = objective_rlm
)
  
rooted_trees <- list()
top_n = 3
for (i in seq_along(objective)) {
  rtree <- root_rtt(
    t = tree,
    tip.dates = tip_dates,
    topx = top_n, 
    ncpu = threads,
    objective = names(objective)[[i]],
    objective_fn = objective[[i]]
  )
  names(rtree) <- paste0(names(objective)[i], "_", 1:top_n)
  index_from = (i-1)*top_n+1
  index_to = i*top_n
  rooted_trees[index_from:index_to] <- rtree
  names(rooted_trees)[index_from:index_to] <- names(rtree)
}

# calculate snps for each rooted tree
snp <- lapply(rooted_trees, function(x) {
  ape::node.depth.edgelength(x)[1:ape::Ntip(tree)]
})

# rescale tip_dates to calendar dates
tip_dates <- as.Date(tip_dates, origin = "1970-01-01")

# recalculate root-to-tip regression on best trees
fit <- lapply(snp, function(x) lm(x~tip_dates))

results <- data.frame(
  r.squared = sapply(fit, function(x) summary(x)$r.squared),
  adj.r.squared = sapply(fit, function(x) summary(x)$r.squared),
  rse = sapply(fit, function(x) summary(x)$sigma),
  ssr = sapply(fit, function(x) sum((summary(x)$residuals)^2)),
  mrca = sapply(fit, function(x) -x$coef[1]/x$coef[2]),
  first = min(tip_dates)[1]
)
results$first <- as.Date(results$first, origin = "1970-01-01")
results$mrca <- as.Date(results$mrca, origin = "1970-01-01")

df <- data.frame()
for (i in seq_along(fit)) {
  new_df <- data.frame(
    name = names(rooted_trees)[i],
    snp = fit[[i]]$model$x,
    date = fit[[i]]$model$tip_dates
  )
  df <- dplyr::bind_rows(df, new_df)
  
}

g <- ggplot(df, aes(date, snp)) + 
  geom_point() + 
  facet_grid(name~.) + 
  geom_smooth(method = "lm")

if (!interactive()) {
  rooted_trees_path <- "rooted_trees.rds"
  rtt_metrics_path <- "rtt_metrics.rds"
  rtt_plots_path <- "rtt_plots.pdf"
  nwk_dir <- "rooted_trees"
} else {
  rooted_trees_path <- paste0(test_dir, "/rooted_trees.rds")
  rtt_metrics_path <- paste0(test_dir, "/rtt_metrics.rds")
  rtt_plots_path <- paste0(test_dir, "/rtt_plots.pdf")
  nwk_dir <- paste0(test_dir, "/rooted_trees")
}

# export rooted tree object
saveRDS(rooted_trees, file = rooted_trees_path)
# export root to tip metrics
saveRDS(results, file = rtt_metrics_path)
# export root to tip regression plots
ggsave(
  filename = rtt_plots_path,
  plot = g,
  width = 10,
  height = 5 * length(rooted_trees),
  limitsize = FALSE
)

if (!dir.exists(nwk_dir)) dir.create(nwk_dir, recursive = TRUE)

for (i in seq_along(rooted_trees)) {
  ape::write.tree(
    rooted_trees[i],
    file = paste0(nwk_dir, "/rooted_tree_rtt_", names(rooted_trees)[i], ".tre")
  )
}

# end logging
if (!interactive()) {
  sink(con)
}
