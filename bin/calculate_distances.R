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
  #####TODO######
  simtree_paths
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
}

simtree_paths <- readLines(simtree_paths)
simtrees <- list()
simtrees$trees <- list()
for (i in simtree_paths) {
  newtrees <- readRDS(i)
  simtrees$trees <- c(simtrees$trees, newtrees$trees)
}

# geographic distance - same country

same_country <- matrix(0, nrow(assemblies), nrow(assemblies))
for (i in unique(assemblies$country)){
  index <- which(assemblies$country == i)
  same_country[index, index] <- 1
}
# mask matrix if analysis is "focused"
if (!is.null(focus_on) & !is.null(focus_by)) {
  same_country <- same_country * maskmat
}
diag(same_country)<-NA
# export data
if(!interactive()) {
  saveRDS(same_country, file = "same_country.rds")
}

# geographic distance - neighbors
neighbors <- matrix(0, nrow(assemblies), nrow(assemblies))
for (i in unique(assemblies$country_iso2c)){
  index1 = which(assemblies$country_iso2c == i)
  data("custom_country_borders")
  borders <- edit_borders(custom_country_borders)
  index2 = which(assemblies$country_iso2c %in% all_neighbors(i, borders = borders))
  
  neighbors[index1, index2] = 1
  neighbors[index2, index1] = 1
}
# mask matrix if analysis is "focused"
if (!is.null(focus_on) & !is.null(focus_by)) {
  neighbors <- neighbors * maskmat
}
diag(neighbors)<-NA
# export data
if (!interactive()) {
  saveRDS(neighbors, file = "neighbors.rds")
}

# geographic distance - same continent

same_continent <- matrix(0, nrow(assemblies), nrow(assemblies))
for (i in unique(assemblies$continent)){
  index = which(assemblies$continent == i)
  same_continent[index, index] = 1
}
# mask matrix if analysis is "focused"
if (!is.null(focus_on) & !is.null(focus_by)) {
  same_continent <- same_continent * maskmat
}
diag(same_continent)<-NA
# export data
if (!interactive()) {
  saveRDS(same_continent, file = "same_continent.rds")
}

# geographic distance - distances in km
geodist <- matrix(NA, nrow(assemblies), nrow(assemblies))
# TODO checking for "lat" and "lon" should be in input validation
if ("lat" %in% names(assemblies) & "lon" %in% names(assemblies)) {
  indices <- 1:nrow(assemblies)
  for (i in 1:nrow(assemblies)) {
    lat1 <- assemblies$lat[i]
    lon1 <- assemblies$lon[i]
    
    # if coordinates are identical set distance to 0.
    index <- which(assemblies$lat == lat1 & assemblies$lon == lon1)
    geodist[i, index] = 0
    geodist[index, i] = 0
    
    index_test <- which(is.na(geodist[i, ]))
    
    if (length(index_test) > 0) {
      s <- assemblies[index_test, which(names(assemblies) %in% c("lat", "lon"))]
      s <- dplyr::distinct(s)
      for (j in 1:nrow(s)) {
        geodist_km <- round(geosphere::distHaversine(
          p1 = c(lon1, lat1),
          p2 = c(s$lon[j], s$lat[j])
        )/1000, 0)
        new <- which(assemblies$lon == s$lon[j] & assemblies$lat == s$lat[j])
        geodist[index, new] <- geodist_km
        geodist[new, index] <- geodist_km
      }
    }
  }
  # mask matrix if analysis is "focused"
  if (!is.null(focus_on) & !is.null(focus_by)) {
    geodist <- geodist * maskmat
  }
  diag(geodist) = NA
}
# export data
if (!interactive()) {
  saveRDS(geodist, file = "geodist.rds")
}

# temporal distance - time difference between collections dates

dates <- unname(lubridate::decimal_date(date_middle(assemblies$collection_date)))
colldist = round(abs(outer(dates, dates, "-")),2)
# mask matrix if analysis is "focused"
if (!is.null(focus_on) & !is.null(focus_by)) {
  colldist <- colldist * maskmat
}
diag(colldist) <- NA
rownames(colldist) <- assemblies$assembly
colnames(colldist) <- assemblies$assembly
# export data
if (!interactive()) {
  saveRDS(colldist, file = "colldist.rds")
}

# temporal distance - cophenetic (patristic) distance between isolates

nsim <- length(simtrees$trees)
ntips <- length(simtrees$trees[[1]]$tip.label)
phylodist_list <- list()

for (i in 1:nsim) {
  phylodist <- ape::cophenetic.phylo(simtrees$trees[[i]])
  index_phylodist <- order(colnames(phylodist))
  phylodist <- phylodist[index_phylodist, index_phylodist]
  # mask matrix if analysis is "focused"
  if (!is.null(focus_on) & !is.null(focus_by)) {
    phylodist <- phylodist * maskmat
  }
  diag(phylodist) <- NA
  
  index_colldist <- unname(sapply(colnames(phylodist), function(x) {
    which(colnames(colldist) == x)
  }))
  colldist_subset <- colldist[index_colldist, index_colldist]

  # mrca definition: distance from the older sample 
  mrca = (phylodist - colldist_subset)/2
  mrca[which(mrca < 0)] = 0
  phylodist_list[[i]] <- mrca
}
# export data
if (!interactive()) {
  saveRDS(phylodist_list, file = "phylodist_list.rds")
}
