#' Plot character state changes over phylogenetic distance
#' 
#' This convenience function integrates a number of other functions in the
#' package. It calculates most recent common ancestors and compares character
#' states between tips of a phylogenetic tree, and then plots the relative
#' frequency of character state changes over phylogenetic distance.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param state_var character; the name of the variable which stores the
#' character state.
#' @param ... additional arguments used with \code{mrca()}.
#' @return a ggplot.
#' @import dplyr
#' @export
plot_statediff_mrca <- function(tree_tbl,
                                state_var,
                                tips = NULL,
                                binwidth = 2,
                                ...) {
  mrca <- mrca(tree_tbl, ...)
  statediff <- statediff(tree_tbl, state_var = state_var)
  df <- dplyr::full_join(mrca, statediff, by = c("tip1", "tip2"))
  if (is.null(tips)) {
    df2 <- df
  } else {
    df2 <- df[which(df$tip1 %in% tips & df$tip2 %in% tips),]
  }
  mrca_max <- max(df2$mrca, na.rm = TRUE)
  internal_breaks <- seq(from = 0, to = mrca_max, by = binwidth)[-1]
  df2$mrca_cat <- cut(df2$mrca, breaks = c(-0.1, internal_breaks, mrca_max))
  df3 <- df2 %>% group_by(mrca_cat) %>% summarise(p = mean(get(state_var)))
  xlab_extra <- switch(attr(mrca, "mrca_method"),
                       first = "(from earlier sample, years)",
                       last = "(from more recent sample, years)",
                       average = "(average, years)")
  g <- ggplot(df3, aes(mrca_cat, p)) + 
    geom_point() +
    xlab(paste0("Most recent common ancestor ", xlab_extra)) +
    ylab(paste0("Relative frequency of ", state_var, " changes")) +
    theme(axis.text.x = element_text(angle = 90))
  g
}
