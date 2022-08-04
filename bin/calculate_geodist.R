rm(list=ls())
library(devtools)
library(dplyr)
library(geosphere)
load_all()

args <- commandArgs(trailingOnly = TRUE)

assemblies <- read.csv(args[1], sep = "\t")
assemblies <- assemblies[order(assemblies$assembly),]

# geographic distance - same country

same_country <- matrix(0, nrow(assemblies), nrow(assemblies))
if ("country" %in% names(assemblies)) {
  for (i in unique(assemblies$country)){
    index = which(assemblies$country == i)
    same_country[index, index] = 1
  }
  diag(same_country)<-NA
}
saveRDS(same_country, file = "same_country.rds")

# geographic distance - neighbors

neighbors <- matrix(0, nrow(assemblies), nrow(assemblies))
if ("iso2c" %in% names(assemblies)) {
  for (i in unique(assemblies$iso2c)){
    index1 = which(assemblies$iso2c == i)
    
    data("custom_country_borders")
    borders <- edit_borders(custom_country_borders)
    index2 = which(assemblies$iso2c %in% all_neighbors(i, borders = borders))
    
    neighbors[index1, index2] = 1
    neighbors[index2, index1] = 1
  }
  diag(neighbors)<-NA
}
saveRDS(neighbors, file = "neighbors.rds")

# geographic distance - same continent

same_continent <- matrix(0, nrow(assemblies), nrow(assemblies))
if ("continent" %in% names(assemblies)) {
  for (i in unique(assemblies$continent)){
    index = which(assemblies$continent == i)
    same_continent[index, index] = 1
  }
  diag(same_continent)<-NA
}
saveRDS(same_continent, file = "same_continent.rds")

# geographic distance - distances in km

geodist <- matrix(NA, nrow(assemblies), nrow(assemblies))
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
  diag(geodist) = NA
}
saveRDS(geodist, file = "geodist.rds")
