#' Plot a phlyogenetic tree using the fan layout
#' 
#' This function plots a phylogenetic tree using the fan layout and optionally
#' adds a number of heatmaps.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param drop_tip character; a vector of tip labels to exclude from the tree.
#' @param open_angle numeric; open angle for the fan layout
#' @param heatmap_var character; a vector or variable names used for creating
#' heatmaps. If \code{NULL} the plot will not contain any heatmaps
#' @param heatmap_colors list; a list of data frames used for color coding
#' heatmaps. See Details for more information. If \code{NULL}, the colors will
#' be generated automatically.
#' @param export logical; should the plot be exported?
#' @details When one or more heatmaps are added to the plot, it is possible to
#' define heatmap colors manually, using the \code{heatmap_colors} argument.
#' By default, the argument is \code{NULL} and the colors are defined
#' automatically. If the argument is not \code{NULL}, the function expects a
#' list of data frames, one data frame for each variable with manually defined
#' colors. The list elements must be named such that the name of the list
#' element matches the variable to which the colors belong.
#' @details. These data frames must contain a unique set of character levels and
#' corresponding color codes. The data frames must contain at least two columns:
#' one for the character levels, and one for the colors. The names of these
#' columns should be the name of the variable, and the \code{"color"},
#' respectively.
#' @examples 
#' \dontrun{
#' # manually define colors for the variable "mlst"
#' mlst_colors <- data.frame(
#'   mlst = c("ST1", "ST2, "Other"),
#'   color = c("red, "green", "grey50")
#' )
#' # plot the phylogenetic tree with manually defined colors
#' plot_tree_fan(
#'   tree_tbl,
#'   heatmap_var = "mlst",
#'   heatmap_colors = "mlst_colors"
#' )
#' }
#' @import ggnewscale
#' @import ggimage
#' @import ggplot2
#' @import ggtree
#' @importFrom qualpalr qualpal
#' @export
plot_tree_fan <- function(tree_tbl,
                          drop_tip = NULL,
                          open_angle = 10,
                          heatmap_var = NULL,
                          heatmap_colors = NULL,
                          heatmap_offset = 0,
                          heatmap_width = 0.1,
                          heatmap_colnames_position = "bottom",
                          heatmap_colnames_angle = 45,
                          heatmap_colnames_offset_x = 0,
                          heatmap_colnames_offset_y = 0,
                          heatmap_font_size = 4,
                          heatmap_hjust = 1,
                          scale = 1,
                          file_name = "tree.pdf",
                          verbose = getOption("verbose")) {
  if (!is.null(file_name) && !grepl("\\.pdf$", file_name)) {
    stop("'filename' must be pdf.")
  }
  if (!is.null(heatmap_var)) {
    if (length(heatmap_var) != length(heatmap_offset)) {
      stop("'heatmap_var' and 'heatmap_offset' must have the same length.")
    }
  }
  tree <-  treeio::as.treedata(tree_tbl)
  tree <- ape::as.phylo(tree)
  if (!is.null(drop_tip)) {
    tree <- ape::drop.tip(tree, tip = drop_tip)
  }
  options(ignore.negative.edge=TRUE)
  p <- ggtree(tree, layout = "fan", open.angle = open_angle)
  if (!is.null(heatmap_var)) {
    idx_x <- which(tree_tbl$label %in% tree$tip.label)
    for (i in 1:length(heatmap_var)) {
      idx_y <- which(names(tree_tbl) == heatmap_var[i])
      hmdf <- as.data.frame(tree_tbl[idx_x, idx_y])
      rownames(hmdf) <- tree_tbl$label[idx_x]
      hmcolors <- "not_set"
      if (!is.null(heatmap_colors)) {
        index <- which(names(heatmap_colors) == heatmap_var[i])
        if (length(index) == 1) {
          if (verbose) {
            message(paste0(
              "Color coding table found for variable ", heatmap_var[i], "."))
          }
          # TODO: Add validation and informative messages around formatting the
          # heatmap color tables
          hmdf_colors <- heatmap_colors[[index]]$color
          names(hmdf_colors) <- heatmap_colors[[index]][[heatmap_var[i]]]
          hmcolors <- "set"
        }
        if (length(index) > 1) {
          stop(paste0(
            "Multiple color coding tables found for variable ", heatmap_var[i])
          )
        }
      }
      if (hmcolors == "not_set") {
        if (verbose) {
          message(paste0(
            "Color coding table not found for variable ",
            heatmap_var[i],
            ". Colors will be assigned automatically."
          ))
        }
        hmdf_colors <- qualpalr::qualpal(
          length(unique(hmdf[[heatmap_var[i]]])), colorspace = "pretty")$hex
        names(hmdf_colors) <- unique(hmdf[[heatmap_var[i]]])
        hmcolors <- "set"
      }
      if (i > 1) {
        p <- p + new_scale_fill()
      }
      p <- gheatmap(
        p,
        hmdf,
        offset = heatmap_offset[i],
        width = heatmap_width,
        colnames_position = heatmap_colnames_position,
        colnames_angle = heatmap_colnames_angle,
        colnames_offset_x = heatmap_colnames_offset_x,
        colnames_offset_y = heatmap_colnames_offset_y,
        font.size = heatmap_font_size,
        hjust = heatmap_hjust
      )+
        scale_fill_manual(
          name = heatmap_var[i],
          values = hmdf_colors,
          breaks = names(hmdf_colors))
    }
  }
  if (!is.null(file_name)) {
    if (!dir.exists(dirname(file_name))) {
      dir.create(dirname(file_name), recursive = TRUE)
    }
    ggsave(
      filename = file_name,
      limitsize = FALSE
    )
    if (verbose) {
      message(paste0("Plot exported to ", file_name, "."))
    }
  } else {
    p
  }
}
