#' Plot a phlyogenetic tree using the fan layout
#' 
#' This function plots a phylogenetic tree using the fan layout and optionally
#' adds a number of heatmaps.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param linewidth numeric; line width for the phylogenetic tree.
#' @param highlight_var character; a variable name used for highlighting tips.
#' @param open_angle numeric; open angle for the fan layout.
#' @param heatmap_var character; a vector or variable names used for creating
#' heatmaps. If \code{NULL} the plot will not contain any heatmaps.
#' @param heatmap_colors list; a list of data frames used for color coding
#' heatmaps. See Details for more information. If \code{NULL}, the colors will
#' be generated automatically.
#' @param heatmap_offset numeric; distance between heatmaps.
#' @param heatmap_width numeric; width of a heatmap band.
#' @param heatmap_colnames_angle numeric; angle for column names.
#' @param heatmap_colnames_font_size numeric; font size for column names.
#' @param heatmap_colnames_hjust numeric; hjust for column names.
#' @param legend_show character; name of heatmap variables to include in legend.
#' By default, all variables are included.
#' @param legend_breaks list; a list of character vectors that tell which
#' levels to show on each heatmap. By default all levels will be shown for
#' each heatmap.
#' @param legend_guides list; a list of lists where the name of each list
#' element is the name of a legend we want to customize and the value is a list
#' of parameters that will be evaluated by guide_legend(). If left empty, the
#' legends will be drawn with defaults.
#' @param legend_box character; arrangement of multiple legends. Can be either
#' \code{"horizontal"} or \code{"vertical"}.
#' @param verbose logical; should verbose messages be printed to the console?
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
#' @import ggtreeExtra
#' @importFrom qualpalr qualpal
#' @export
plot_tree_fan <- function(tree_tbl,
                          linewidth = 0.5,
                          highlight_var = NULL,
                          open_angle = 10,
                          heatmap_var = NULL,
                          heatmap_colors = NULL,
                          heatmap_offset = 0.15,
                          heatmap_width = 5,
                          heatmap_colnames_angle = 45,
                          heatmap_colnames_font_size = 4,
                          heatmap_colnames_hjust = 1,
                          legend_show = NA,
                          legend_breaks = list(),
                          legend_guides = list(),
                          legend_box = "horizontal",
                          file_name = NULL,
                          verbose = getOption("verbose")) {
  legend_box <- match.arg(legend_box, choices = c("horizontal", "vertical"))
  if (!is.null(file_name) && !grepl("\\.pdf$", file_name)) {
    stop("'filename' must be pdf.")
  }
  tree_db <-  treeio::as.treedata(tree_tbl)
  tree <- ape::as.phylo(tree_db)
  options(ignore.negative.edge=TRUE)
  if (is.null(highlight_var)) {
    p <- ggtree(
      tree_db,
      size = linewidth,
      layout = "fan",
      open.angle = open_angle
    )
  } else {
    p <- ggtree(
      tree_db,
      aes(color = get(highlight_var)),
      size = linewidth,
      layout = "fan",
      open.angle = open_angle
    ) + 
    labs(
      color = highlight_var
    )
    if (!highlight_var %in% legend_show){
      p <- p + guides(
        color = "none"
      )
    }
  }
  if (!is.null(heatmap_var)) {
    idx_x <- which(tree_tbl$label %in% tree$tip.label)
    for (i in 1:length(heatmap_var)) {
      # define heatmap data frame
      idx_y <- which(names(tree_tbl) == heatmap_var[i])
      hmdf <- as.data.frame(tree_tbl[idx_x, idx_y])
      row.names(hmdf) <- tree_tbl$label[idx_x]
      hmdf <- data.frame(id = row.names(hmdf), group = hmdf[[heatmap_var[i]]])
      # define heatmap colors
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
      # add heatmaps
      p <- p + new_scale_fill()
      axis_text <- heatmap_var[i]
      if (length(legend_show) == 1 && is.na(legend_show)) {
        show_legend <- NA
      } else {
        show_legend <- ifelse(heatmap_var[i] %in% legend_show, NA, FALSE)
      }
      p <- p + geom_fruit(
        data=hmdf,
        geom=geom_tile,
        mapping=aes(y=id, fill=group),
        width=heatmap_width,
        offset=heatmap_offset,
        axis.params = list(axis = "x",
                           text = axis_text,
                           text.size = heatmap_colnames_font_size,
                           text.angle = heatmap_colnames_angle,
                           hjust = heatmap_colnames_hjust,
                           line.size = 0,
                           line.color = "#FFFFFF"),
        show.legend = show_legend
      )
      
      
      if (heatmap_var[i] %in% names(legend_breaks)) {
        breaks <- legend_breaks[[heatmap_var[i]]]
      } else {
        breaks <- names(hmdf_colors)
      }
      
      if (heatmap_var[i] %in% names(legend_guides)) {
        p <- p + scale_fill_manual(
          name = heatmap_var[i],
          values = hmdf_colors,
          breaks = breaks,
          na.translate = FALSE,
          guide = rlang::exec(guide_legend, !!!legend_guides[[heatmap_var[i]]])
        )
      } else {
        p <- p + scale_fill_manual(
          name = heatmap_var[i],
          values = hmdf_colors,
          breaks = breaks,
          na.translate = FALSE
        )
      }
      
      
    }
  }
  p <- p + theme(legend.box = legend_box)
  p
}
