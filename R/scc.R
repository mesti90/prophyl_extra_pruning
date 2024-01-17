#' Count state changes between two states
#' @param target character; the target ancestral state to plot. Must match a
#' column name in \code{tree_tbl}.
#' @param from character; origin
#' @param to character; destination
#' @param tree_tbl tibble; dated tree with predicted ancestral states. This
#' tibble is automatically returned by the nextlow pipeline.
#' @examples 
#' \dontrun {
#' scc(from = "KL3", to = "KL2", target = "k_serotype", tree_tbl = tree_tbl)
#' }
#' @export
scc <- function(from, to, target, tree_tbl) {
  target <- match.arg(target, names(tree_tbl))
  cats <- unique(unlist(strsplit(tree_tbl[[target]], split = "\\|")))
  from <- match.arg(from, cats)
  to <- match.arg(to, cats)
  match_from <- mapply(function(x, y) {
    hits <- strsplit(x, split = "\\|")[[1]]
    y %in% hits
  }, tree_tbl[[paste0(target, "_from")]], from)
  match_to <- mapply(function(x, y) {
    hits <- strsplit(x, split = "\\|")[[1]]
    y %in% hits
  }, tree_tbl[[target]], to)
  stree <- tree_tbl[which(match_from & match_to),]
  scprob <- unname(unlist(mapply(function(x, y) {
    scp(x, y, from_which = from, to_which = to)$prob
  }, stree[[paste0(target, "_from")]], stree[[target]])))
  out <- list(tree_tbl = stree, prob = scprob, count = sum(scprob))
  return(out)
}
