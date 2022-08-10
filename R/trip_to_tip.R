#' Trip to tip analysis
#' 
#' Trip to tip analysis counts the number of ancestral state changes between a
#' common ancestor and each tip which descends from this common ancestor. The 
#' function currently counts the number of state changes from root. 
#' @param tree_tbl tibble
#' @import ggplot2
#' @examples 
#' \dontrun{
#' res <- trip_to_tip(tree_tbl, "k_serotype", min_sc_prob = 0.5)
#' }
#' @export
trip_to_tip <- function(tree_tbl, target, focus = NULL, min_sc_prob = 0) {
  index_root <- which(tree_tbl$parent == tree_tbl$node)
  if (length(index_root) == 0) {
    stop("Root node seems to be missing. Check.")
  }
  if (length(index_root) > 1) {
    stop("Multiple root nodes detected. Check.")
  }
  tree_tbl$collection_day <- lubridate::decimal_date(tree_tbl$collection_day)
  root_day <- tree_tbl$collection_day[index_root]
  tree_tbl$collection_day <- tree_tbl$collection_day-root_day
  fcol <- paste0(target, "_sc_prob")
  if (fcol %in% names(tree_tbl) == FALSE) {
    stop(paste0("Variable ", target, "_sc_prob not found. Check."))
  }
  pcol <- paste0(target, "_sc")
  if (pcol %in% names(tree_tbl) == FALSE) {
    stop(paste0("Variable ", target, "_sc not found. Check."))
  }
  tree_tbl$position <- ifelse(grepl("^Node_", tree_tbl$label), "inode", "tip")
  if (!is.null(focus)) {
    if (mean(focus %in% tree_tbl[[target]]) < 1) {
      miss <- focus[which(focus %in% tree_tbl[[target]] == FALSE)]
      miss <- paste(miss, collapse = ", ")
      stop(paste0("focus = ", miss, " not found within ", target, ". Check."))
    } else {
      tree_tbl[[target]] <- ifelse(tree_tbl[[target]] %in% focus,
                                   tree_tbl[[target]], 
                                   "Other")
    }
  }
  index <- which(tree_tbl[[fcol]] >= min_sc_prob)
  tree_tbl2 <- tree_tbl[index, ]
  root_day <- tree_tbl2$collection_day[tree_tbl2$label == "Node_1"]
  # fit linear model for tips
  fit_tip_x <- tree_tbl2$collection_day[which(tree_tbl2$position == "tip")]
  fit_tip_y <- tree_tbl2[[pcol]][which(tree_tbl2$position == "tip")]
  fit_tip <- lm(fit_tip_y ~ fit_tip_x + 0)
  tip_msg <- paste0(
    "y = ", round(fit_tip$coeff[1], 3), "x \n",
    "R^2 = ", round(summary(fit_tip)$r.squared, 3)
  )
  # fit linear model for internal nodes
  fit_inode_x <- tree_tbl2$collection_day[which(tree_tbl2$position == "inode")]
  fit_inode_y <- tree_tbl2[[pcol]][which(tree_tbl2$position == "inode")]
  fit_inode <- lm(fit_inode_y ~ fit_inode_x + 0)
  inode_msg <- paste0(
    "y = ", round(fit_inode$coeff[1], 3), "x \n",
    "R^2 = ", round(summary(fit_inode)$r.squared, 3)
  )
  # plot trip-to-tip
  xmax <- max(tree_tbl2$collection_day, na.rm = TRUE)
  ymax <- max(tree_tbl2[[pcol]], na.rm = TRUE)
  cols <- qualpalr::qualpal(2, colorspace = "pretty")$hex
  if (is.null(focus)) {
    g <- ggplot(tree_tbl2, aes(collection_day, get(pcol)))+
      geom_point(aes(col = position,
                     alpha = position),
                 size = 4)+
      geom_smooth(method = lm,
                  formula = "y ~ x + 0",
                  se = FALSE,
                  fullrange = TRUE,
                  aes(col = position),
                  linetype = "dotted")+
      scale_color_manual(values = cols)+
      scale_alpha_manual(breaks = c("tip", "inode"), values = c(1, 0.5))+
      xlim(0, xmax)+
      xlab("Time from root (years)")+
      ylab(paste0(target, " changes (count)"))+
      annotate("text",
               x = 0.1*xmax,
               y = 0.85*ymax,
               label = tip_msg,
               hjust = 0,
               col = cols[1])+
      annotate("text",
               x = 0.1*xmax,
               y = 0.65*ymax,
               label = inode_msg,
               hjust = 0,
               col = cols[2])
  } else {
    g <- ggplot(tree_tbl2, aes(collection_day, get(pcol)))+
      geom_point(aes(col = position,
                     alpha = position,
                     shape = get(target),
                     size = get(target)))+
      geom_smooth(method = lm,
                  formula = "y ~ x + 0",
                  se = FALSE,
                  fullrange = TRUE,
                  aes(col = position),
                  linetype = "dotted")+
      scale_color_manual(values = cols)+
      scale_alpha_manual(breaks = c("tip", "inode"), values = c(1, 0.5))+
      scale_size_manual(breaks = c(focus, "Other"),
                        values = c(rep(5, times = length(focus)), 2))+
      xlim(0, xmax)+
      xlab("Time from root (years)")+
      ylab(paste0(target, " changes (count)"))+
      annotate("text",
               x = 0.1*xmax,
               y = 0.85*ymax,
               label = tip_msg,
               hjust = 0,
               col = cols[1])+
      annotate("text",
               x = 0.1*xmax,
               y = 0.65*ymax,
               label = inode_msg,
               hjust = 0,
               col = cols[2])
  }
  
  return(list(fit_tip = fit_tip, fit_inode = fit_inode, plot = g))
}
