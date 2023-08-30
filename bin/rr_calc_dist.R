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

library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "Path to the assemblies file."
  ),
  make_option(
    c("-s", "--simtrees"),
    type = "character",
    help = "Path to a text file containing simtree paths."
  ),
  make_option(
    c("-b", "--focus_by"),
    type = "character",
    help = "Variable to focus on."
  ),
  make_option(
    c("-o", "--focus_on"),
    type = "character",
    help = "Value of the variable to focus on."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/prophyl",
    assemblies = "assemblies.tsv",
    simtrees = "simtree_paths.txt",
    focus_by = "none",
    focus_on = "none"
  )
}

library(devtools)
library(dplyr)
library(geosphere)
library(lubridate)

load_all(args$project_dir)

assemblies <- read.csv(args$assemblies, sep = "\t")
assemblies <- assemblies[order(assemblies$assembly),]

focus_by <- args$focus_by
focus_on <- args$focus_on

if (focus_by == "none") focus_by <- NULL
if (focus_on == "none") focus_on <- NULL

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

simtree_paths <- readLines(args$simtrees)
simtrees <- list()
simtrees$trees <- list()
for (i in simtree_paths) {
  newtrees <- readRDS(i)
  simtrees$trees <- c(simtrees$trees, newtrees$trees)
}

# geographic distance - same city
same_city <- varid_matrix(
  df = assemblies,
  var = "city",
  focus_by = focus_by,
  focus_on = focus_on
)

# export data
saveRDS(same_city, file = "same_city.rds")

# geographic distance - same country
same_country <- varid_matrix(
  df = assemblies,
  var = "country",
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(same_country, file = "same_country.rds")

# geographic distance - neighbors
neighbors <- neighbors_matrix(
  df = assemblies,
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(neighbors, file = "neighbors.rds")

# geographic distance - same continent
same_continent <- varid_matrix(
  df = assemblies,
  var = "continent",
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(same_continent, file = "same_continent.rds")

# geographic distance - distances in km
geodist <- geodist_matrix(
  df = assemblies,
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(geodist, file = "geodist.rds")

# temporal distance - time difference between collections dates
colldist <- colldist_matrix(
  df = assemblies,
  focus_by = focus_by,
  focus_on = focus_on,
  estimate_dates = "middle"
)
# export data
saveRDS(colldist, file = "colldist.rds")

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
saveRDS(phylodist_list, file = "phylodist_list.rds")
