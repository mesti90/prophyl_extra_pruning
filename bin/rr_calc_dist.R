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
    c("-t", "--tree"),
    type = "character",
    help = "A dated tree in rds format."
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
    tree = "dated_tree.rds",
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

# geographic distance - same city
same_city <- varid_matrix(
  df = assemblies,
  id_var = "assembly",
  var = "city",
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(same_city, file = "same_city.rds")

# geographic distance - same country
same_country <- varid_matrix(
  df = assemblies,
  id_var = "assembly",
  var = "country",
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(same_country, file = "same_country.rds")

# geographic distance - neighbors

data("custom_country_borders")
updated_borders <- edit_borders(custom_country_borders)

neighbors <- neighbors_matrix(
  df = assemblies,
  id_var = "assembly",
  iso2c_var = "country_iso2c",
  country_borders = updated_borders,
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(neighbors, file = "neighbors.rds")

# geographic distance - same continent
same_continent <- varid_matrix(
  df = assemblies,
  id_var = "assembly",
  var = "continent",
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(same_continent, file = "same_continent.rds")

# geographic distance - distances in km
geodist <- geodist_matrix(
  df = assemblies,
  id_var = "assembly",
  lat_var = "lat",
  lon_var = "lon",
  focus_by = focus_by,
  focus_on = focus_on
)
# export data
saveRDS(geodist, file = "geodist.rds")

# temporal distance - time difference between collection dates

# set random seed for date estimation
set.seed(0)

colldist <- tempdist_matrix(
  df = assemblies,
  id_var = "assembly",
  date_var = "collection_date",
  focus_by = focus_by,
  focus_on = focus_on,
  estimate_dates = "runif"
)
# export data
saveRDS(colldist, file = "colldist.rds")

# temporal distance - most recent common ancestors between isolates

# for the big tree

tree <- readRDS(args$tree)

phylodist <- phylodist_matrix(
    tree = tree,
    df = assemblies,
    id_var = "assembly",
    focus_by = focus_by,
    focus_on = focus_on
)

# shrink colldist to phylodist_subset, align rows and columns
colldist_subset <- shrink_matrix(colldist, phylodist)
  
# calculate mrca
mrca <- mrca_matrix(
  phylodist,
  colldist_subset,
  force_nonnegative = TRUE
)

# export data
saveRDS(mrca, file = "phylodist.rds")

# for the subsampled trees

simtree_paths <- readLines(args$simtrees)
simtree_names <- gsub("\\.rds","",basename(simtree_paths))

simtrees <- list()
simtrees$trees <- list()
for (i in simtree_paths) {
  newtrees <- readRDS(i)
  simtrees$trees <- c(simtrees$trees, newtrees$trees)
}

phylodist_list <- list()
for (i in seq_along(simtrees$trees)) {
  phylodist_subset <- phylodist_matrix(
    tree = simtrees$trees[[i]],
    df = assemblies,
    id_var = "assembly",
    focus_by = focus_by,
    focus_on = focus_on
  )
  
  # shrink colldist to phylodist_subset, align rows and columns
  colldist_subset <- shrink_matrix(colldist, phylodist_subset)
  
  # calculate mrca
  mrca <- mrca_matrix(
    phylodist_subset,
    colldist_subset,
    force_nonnegative = TRUE
  )
  phylodist_list[[i]] <- mrca
}
names(phylodist_list) <- simtree_names
# export data
saveRDS(phylodist_list, file = "phylodist_list.rds")
