rm(list = ls())
library(igraph)
library(leaflet)
library(leaflet.minicharts)

args <- commandArgs(trailingOnly = TRUE)
jobname <- args[1]
varname <- args[2]

tree <- ape::read.tree(
  paste0("./jobs/", jobname, "/output/treedater/treedater_tree.nwk"))

load(paste0("./jobs/", jobname, "/treemeta.rda"))
aci <- dplyr::rename(aci, label = assembly)

set.seed(0)
o <- unname(sapply(tree$tip.label, function(x) which(aci$label == x)))
prop <- aci$country[o]
names(prop) <- aci$label[o]
tree$edge.length <- ifelse(
  tree$edge.length > 0, tree$edge.length, abs(jitter(tree$edge.length)))
ans <- ape::ace(prop, tree, type = "d")

tree_tbl <- tibble::as_tibble(tree)
index_na <- which(is.na(tree_tbl$label))

aux <- dplyr::left_join(tree_tbl[-index_na, ], aci[,c(
  "label",
  varname,
  "collection_year",
  "collection_day")], by = "label")

tree_tbl <- dplyr::bind_rows(aux, tree_tbl[index_na,])

# ancestral state is the state with the highest likelihood
i <- apply(ans$lik.anc, 1, function(x) {
  index <- which(x == max(x))
  if (length(index == 1) & x[index] >= 0.5) {
    return(index)
  } else {
    return(NA)
  }
})
tree_tbl$prop <- c(prop, sort(unique(prop))[i])
tree_tbl$prop <- ifelse(is.na(tree_tbl$prop), "??", tree_tbl$prop)

n <- length(unique(tree_tbl$prop))-1
adjmat <- matrix(0, nrow = n, ncol = n)
rownames(adjmat) <- sort(unique(tree_tbl$prop))[-1]
colnames(adjmat) <- sort(unique(tree_tbl$prop))[-1]

for (i in 1:nrow(tree_tbl)) {
  prop_child <- tree_tbl$prop[i]
  prop_parent <- tree_tbl$prop[which(tree_tbl$node == tree_tbl$parent[i])]
  if (prop_child == prop_parent) next() else {
    index_x <- which(rownames(adjmat) == prop_parent)
    index_y <- which(colnames(adjmat) == prop_child)
    adjmat[index_x, index_y] <- adjmat[index_x, index_y]+1
  }
}

#make more general
load("country_gps.rda")
country_gps$country <- unique(aci$country)

#set.seed(0)
#adjmat <- matrix(0, nrow = 31, ncol = 31)
#rownames(adjmat) <- sort(unique(aci[[varname]]))
#colnames(adjmat) <- sort(unique(aci[[varname]]))

#for (i in 1:nrow(adjmat)) {
#  for (j in 1:ncol(adjmat)){
#    adjmat[i,j] <- ifelse(i == j, 0, sample(c(rep(0, times = 19), 1), 1))
#  }
#}

edgelist <- reshape2::melt(adjmat)
names(edgelist) <- c("from", "to", "freq")

#test <- igraph::graph_from_adjacency_matrix(adjmat)
#test2 <- get.data.frame(test)
edgelist$from_lat <- sapply(edgelist$from, function(x) {
  country_gps$lat[which(country_gps$country == x)]
})
edgelist$from_lng <- sapply(edgelist$from, function(x) {
  country_gps$lon[which(country_gps$country == x)]
})
edgelist$to_lat <- sapply(edgelist$to, function(x) {
  country_gps$lat[which(country_gps$country == x)]
})
edgelist$to_lng <- sapply(edgelist$to, function(x) {
  country_gps$lon[which(country_gps$country == x)]
})
#edgelist$date <- 0

#edgelist2 <- edgelist
#edgelist2$date <- 1
#edgelist2$freq <- sample(edgelist2$freq, length(edgelist2$freq))

#edgelist <- rbind(edgelist, edgelist2)

#azert nem megy, mert a test2-ben benne kellene lenne annak is, amikor nincs edge?
#0 thickness. minden datum allashoz tartozzon egy graf

map <- leaflet() %>% 
  addTiles() %>%
  addMarkers(lat = country_gps$lat,
             lng = country_gps$lon,
             popup = country_gps$country,
             clusterOptions=markerClusterOptions()) %>%
  addFlows(lat0 = edgelist$from_lat,
           lng0 = edgelist$from_lng,
           lat1 = edgelist$to_lat,
           lng1 = edgelist$to_lng,
           flow = edgelist$freq,
           #time = edgelist$date,
           minThickness = min(adjmat),
           maxThickness = max(adjmat))
library(htmlwidgets)
saveWidget(map, file="map.html")
