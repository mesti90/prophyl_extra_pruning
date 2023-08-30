# This script calculates various matrices for risk calculations. These matrices
# are symmetrical matrices and reflect comparisons between isolates. E.g. the
# matrix "same_country" describes whether two isolates were collected from the
# same country; "geodist" contains geographical distances between two isolates,
# "colldist" contains differences in collection dates, etc. It is possible to 
# run the script without any "focus" in this case all isolates on the tree will
# be included in distance comparisons. Alternatively, it is possible to select a
# focus group. In this case, each isolate can be categorised as "in-focus" or
# "not-in-focus". Comparisons between two "in-focus" isolates will be kept also
# comparisons between an "in-focus" and a "not-in-focus" isolate, but
# comparisons between two "not-in-focus" isolates will be eliminated from the
# analysis. In practice, these elements will be masked with NAs in all matrices.

rm(list=ls())

library(devtools)
library(dplyr)
library(geosphere)
library(lubridate)

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  project_dir <- args[1]
  assemblies_path <- args[2]
  simtree_paths <- args[3]
  focus_by <- args[4]
  focus_on <- args[5]
} else {
  project_dir <- "~/Methods/prophyl"
  test_dir <- "~/Methods/prophyl-tests/test-calculate_distances"
  assemblies_path <- paste0(test_dir, "/assemblies.tsv")
  simtree_paths <- paste0(test_dir, "/test.txt")
  focus_by <- "continent"
  focus_on <- "europe"
}


load_all(project_dir)

assemblies <- read.csv(assemblies_path, sep = "\t")
assemblies <- assemblies[order(assemblies$assembly),]

# TODO maybe this should be checked at input validation as well?
# TODO maybe only character variables should be allowed?
# validate input
if (!is.null(focus_by)) {
  focus_by <- match.arg(focus_by, choices = names(assemblies))
}
if (!is.null(focus_on)) {
  focus_on <- match.arg(focus_on, choices = unique(assemblies[[focus_by]]))
}
if (!is.null(focus_on) & !is.null(focus_by)) {
  maskmat <- mask_matrix(assemblies, focus_by = focus_by, focus_on = focus_on)
  rownames(maskmat) <- assemblies$assembly
  colnames(maskmat) <- assemblies$assembly
}

simtree_paths <- readLines(simtree_paths)
simtrees <- list()
simtrees$trees <- list()
for (i in simtree_paths) {
  newtrees <- readRDS(i)
  simtrees$trees <- c(simtrees$trees, newtrees$trees)
}

# geographic distance - same country
same_country <- varid_matrix(
  df = assemblies,
  var = "country",
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
if(!interactive()) {
  saveRDS(same_country, file = "same_country.rds")
}

# geographic distance - neighbors
neighbors <- neighbors_matrix(
  df = assemblies,
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
if (!interactive()) {
  saveRDS(neighbors, file = "neighbors.rds")
}

# geographic distance - same continent
same_continent <- varid_matrix(
  df = assemblies,
  var = "continent",
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
if (!interactive()) {
  saveRDS(same_continent, file = "same_continent.rds")
}

# geographic distance - distances in km
geodist <- geodist_matrix(
  df = assemblies,
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
if (!interactive()) {
  saveRDS(geodist, file = "geodist.rds")
}

# temporal distance - time difference between collections dates
colldist <- colldist_matrix(
  df = assemblies,
  focus_by = focus_by,
  focus_on = focus_on,
  estimate_dates = "middle"
)
# export data
if (!interactive()) {
  saveRDS(colldist, file = "colldist.rds")
}

# temporal distance - most recent common ancestors between isolates

nsim <- length(simtrees$trees)
phylodist_list <- list()
for (i in 1:nsim) {
  phylodist_subset <- phylodist_matrix(
    tree = simtrees$trees[[i]],
    df = assemblies,
    focus_by = focus_by,
    focus_on = focus_on
  )
  # subset colldist to relevant rows and columns
  index <- unname(sapply(colnames(phylodist_subset), function(x) {
    which(colnames(colldist) == x)
  }))
  colldist_subset <- colldist[index, index]
  # calculate mrca
  mrca <- mrca_matrix(
    phylodist_subset,
    colldist_subset,
    force_nonnegative = TRUE
  )
  phylodist_list[[i]] <- mrca
}
# export data
if (!interactive()) {
  saveRDS(phylodist_list, file = "phylodist_list.rds")
}
