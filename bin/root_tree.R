# This script takes a phylogenetic tree and attempts to root it using three
# aproaches:
#
# 1. Midpoint rooting
# 2. Minimum Ancestor Deviation (MAD)
# 3. Root-to-tip regression
#
# Root-to-tip regression uses a custom function that was built on the
# non-exported .multi.rtt() function from the treedater package which in turn
# was built on the rtt() function from the ape package. For more details, check
# the function documentation ?root_rtt().

library(devtools)
library(optparse)
rm(list = ls())

args_list <- list(
 make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
 make_option(
   c("-t", "--tree"),
   type = "character",
   help = "Path to tree file."
 ),
 make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "Path to assemblies file."
  ),
  make_option(
    c("-c", "--threads"),
    type = "integer",
    help = "Number of threads to use."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/prophyl",
    tree = "treeshrink.tre",
    assemblies = "assemblies.tsv",
    threads = 10
  )
}

load_all(args$project_dir)

# create log file and start logging
if (!interactive()) {
  con <- file("log.txt")
  sink(con, split = TRUE)
}

# read tree
tree <- ape::read.tree(args$tree)

# if tree is rooted, unroot
if (ape::is.rooted(tree)) {
  tree <- ape::unroot(tree)
}

rooted_trees <- list()

# OPTION 1: MIDPOINT ROOTING

# root tree
rooted_trees[["midpoint"]] <- phytools::midpoint.root(tree)

# Note: OPTION 2  is removed until Issue #73 is fixed.
# # OPTION 2: MINIMUM ANCESTOR DEVIATION (MAD)

# # root tree
# mad <- root_mad(
#   tree,
#   output_mode = "full",
#   cache = TRUE,
#   threads = args$threads,
#   verbose = TRUE
# )

# # add tips that were collapsed

# mad_tree <- mad[[3]]
# collapsed_tips <- mad[[7]]

# for (i in 1:nrow(collapsed_tips)) {
#   mad_tree <- TreeTools::AddTip(
#     mad_tree,
#     where = collapsed_tips$keep[i],
#     label = collapsed_tips$drop[i],
#     edgeLength = 0
#   )
# }

# rooted_trees[["mad"]] <- mad_tree

# OPTION 3: ROOT-TO-TIP REGRESSION

# read assemblies
assemblies <- read.csv(args$assemblies, sep = "\t")

# collect tip dates in the same order as tree$tip.label
tip_dates <- unname(sapply(tree$tip.label, function(x) {
  index <- which(assemblies$assembly == x)
  assemblies$collection_day[index]
}))

# convert tip dates to numeric for root_tree()
tip_dates <- as.numeric(as.Date(tip_dates))

# TODO: look for better objectives
objective_rlm_slope <- function(x,y) MASS::rlm(y ~ x)$coef[2]
objective_rlm_rms <- function(x,y) -summary(MASS::rlm(y ~ x))$sigma^2

# Remove custom objectives for now
# objective <- list(
#   "correlation" = NULL,
#   "rsquared" = NULL,
#   "rms" = NULL,
#   "rlm_slope" = objective_rlm_slope,
#   "rlm_rms" = objective_rlm_rms
# )

objective <- list(
  "correlation" = NULL,
  "rsquared" = NULL,
  "rms" = NULL
)

# return the top_n trees for each objective
top_n = 3

for (i in seq_along(objective)) {
  rtree <- root_rtt(
    t = tree,
    tip.dates = tip_dates,
    topx = top_n, 
    ncpu = args$threads,
    objective = names(objective)[[i]],
    objective_fn = objective[[i]]
  )
  names(rtree) <- paste0("rtt_", names(objective)[i], "_", 1:top_n)
  index_from = (i-1)*top_n + 2
  index_to = i*top_n + 1
  rooted_trees[index_from:index_to] <- rtree
  names(rooted_trees)[index_from:index_to] <- names(rtree)
}

# CALCULATE ROOT TO TIP METRICS FOR EACH ROOTED TREE

# calculate snps for each rooted tree
snp <- lapply(rooted_trees, function(x) {
  ape::node.depth.edgelength(x)[1:ape::Ntip(tree)]
})

# rescale tip_dates to calendar dates
tip_dates <- as.Date(tip_dates, origin = "1970-01-01")

# recalculate root-to-tip regression using calendar dates
fit <- lapply(snp, function(x) lm(x~tip_dates))

# calculate metrics for each fit
results <- data.frame(
  r.squared = sapply(fit, function(x) summary(x)$r.squared),
  adj.r.squared = sapply(fit, function(x) summary(x)$adj.r.squared),
  rse = sapply(fit, function(x) summary(x)$sigma),
  ssr = sapply(fit, function(x) sum((summary(x)$residuals)^2)),
  mrca = sapply(fit, function(x) -x$coef[1]/x$coef[2]),
  first = min(tip_dates, na.rm = TRUE)[1]
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

# export rooted tree object
saveRDS(rooted_trees, file = "rooted_trees.rds")
# export root to tip metrics
saveRDS(results, file = "rtt_metrics.rds")
# export root to tip regression plots
ggsave(
  filename = "rtt_plots.pdf",
  plot = g,
  width = 10,
  height = 5 * length(rooted_trees),
  limitsize = FALSE
)

# if (!dir.exists("rooted_trees")) dir.create("rooted_trees")

# for (i in seq_along(rooted_trees)) {
#   ape::write.tree(
#     rooted_trees[i],
#     file = paste0("rooted_trees/rooted_tree_", names(rooted_trees)[i], ".tre")
#   )
# }

# consistency checks

# check that all trees have the same number of tips
ntips <- sapply(rooted_trees, function(x) length(x$tip.label))
testthat::expect_equal(length(unique(ntips)), 1)

# end logging
if (!interactive()) {
  sink(con)
}
