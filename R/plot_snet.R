#' Plot ancestral state change network
#' 
#' This function prepares a plot from an adjacency matrix which represents
#' directed ancestral state network along a phylogenetic tree.
#' @param snet matrix; and adjacency matrix which contains frequencies
#' of ancestral state changes.
#' @param min_freq numeric; the lowest frequency of observed state changes 
#' which should be included on the plot.
#' @param cache_file character; the name of the cache file without the file
#' extension. If \code{NULL}, results are not cached.
#' @import igraph 
#' @examples 
#' \dontrun{
#' plot_snet(snet, min_freq = 3)
#' }
#' @export
plot_snet <- function(snet, min_freq = 0.5) {
  for (i in 1:nrow(snet)){
    for (j in 1:ncol(snet)) {
      snet[i,j] <- ifelse(snet[i,j] >= min_freq, snet[i,j], 0)
    }
  }
  net <- igraph::graph_from_adjacency_matrix(snet)
  net <- igraph::simplify(net, remove.multiple = T, remove.loops = T)
  deg <- igraph::degree(net, mode="all")
  V(net)$size <- 5*log(deg)+2
  l <- igraph::layout_in_circle(net)
  plot(net, edge.arrow.size=.4, layout = l)
}
