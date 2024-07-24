#' Plot a phlyogenetic tree using the fan layout
#' 
#' This function plots a phylogenetic tree using the rectangular layout and
#' optionally adds a number of heatmaps.
#' @param tree_tbl tibble; the phylogenetic tree in tibble format.
#' @param show_tiplab logical; show tip labels?
#' @param tiplab_size logical; text size for tiplabs.
#' @param mrsd tibble; the date of the most recent tip.
#' @param linewidth numeric; line width for the phylogenetic tree.
#' @param highlight_var character; a variable name used for highlighting tips.
#' @param heatmap_var character; a vector or variable names used for creating
#' heatmaps. If \code{NULL} the plot will not contain any heatmaps.
#' @param heatmap_colors list; a list of data frames used for color coding
#' heatmaps. See Details for more information. If \code{NULL}, the colors will
#' be generated automatically.
#' @param heatmap_text boolean; should heatmap values be printed in each heatmap
#' tile?
#' @param heatmap_offset numeric; distance between heatmaps.
#' @param heatmap_width numeric; width of a heatmap band.
#' @param heatmap_colnames_angle numeric; angle for column names.
#' @param heatmap_colnames_font_size numeric; font size for column names.
#' @param heatmap_colnames_hjust numeric; hjust for column names.
#' @param heatmap_colnames_vjust numeric; vjust for column names, down, from
#' bottom.
#' @param legend_show character; name of heatmap variables to include in legend.
#' By default, all variables are included.
#' @param legend_breaks list; a list of character vectors that tell which
#' levels to show on each heatmap. By default all levels will be shown for
#' each heatmap.
#' @param legend_guides list; a list of lists where the name of each list
#' element is the name of a legend we want to customize and the value is a list
#' of parameters that will be evaluated by guide_legend(). If left empty, the
#' legends will be drawn with defaults.
#' @param legend_position character; position of the legends compared to the
#' plots. See `?theme()` for accepted values.
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
#'   mlst = c("ST1", "ST2", "Other"),
#'   color = c("red", "green", "grey50")
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
#' @import patchwork
#' @importFrom qualpalr qualpal
#' @importFrom dplyr %>%
#' @export
plot_tree_long <- function(
  tree_tbl,
  show_tiplab = TRUE,
  tiplab_size = 1,
  mrsd = NULL,
  linewidth = 0.5,
  highlight_var = NULL,
  heatmap_var = NULL,
  heatmap_colors = NULL,
  heatmap_text = TRUE,
  heatmap_offset = 0,
  heatmap_width = 5,
  heatmap_colnames_angle = 0,
  heatmap_colnames_font_size = 4,
  heatmap_colnames_hjust = 1,
  heatmap_colnames_vjust = 0,
  legend_show = NA,
  legend_breaks = list(),
  legend_guides = list(),
  legend_position = "right",
  legend_box = "horizontal",
  verbose = getOption("verbose")
  ) {
  legend_box <- match.arg(legend_box, choices = c("horizontal", "vertical"))
  
  # input validation for mrsd
  if (!is.null(mrsd)) {
    if (length(mrsd) > 1) {
      stop("Argument 'mrsd' must be a single value.")
    } else if(is.na(mrsd)) {
      stop("Argument 'mrsd' cannot be NA.")
    } else if (class(mrsd) != "Date") {
      stop("Argument 'mrsd' must be a 'Date'.")
    }
  }
  
  tree_db <-  treeio::as.treedata(tree_tbl)
  tree <- ape::as.phylo(tree_db)
  options(ignore.negative.edge=TRUE)
  
  # plot base tree, with/without highlighting, with/without mrsd
  if (is.null(highlight_var)) {
    if (is.null(mrsd)) {
      p <- ggtree(tree_db, size = linewidth)
    } else {
      p <- ggtree(tree_db, size = linewidth, mrsd = mrsd) + theme_tree2()
    }
  } else {
    if (is.null(mrsd)) {
      p <- ggtree(tree_db, aes(color = get(highlight_var)), size = linewidth) + 
        labs(
          color = highlight_var
        )
    } else {
      p <- ggtree(
        tree_db, aes(color = get(highlight_var)), size = linewidth, mrsd = mrsd) + 
        labs(
          color = highlight_var
        ) + theme_tree2()
    }
    if (!highlight_var %in% legend_show){
      p <- p + guides(
        color = "none"
      )
    }
  }
  
  # add tiplab
  if (show_tiplab) {
    p <- p + geom_tiplab(
      align = FALSE, geom = "text", size = tiplab_size, hjust = 0)
  }
  
  # add heatmaps
  hm_list <- list()
  if (!is.null(heatmap_var)) {

    # tree_tbl index for tips
    idx_x <- which(tree_tbl$label %in% tree$tip.label)
    
    # from aplot:::get_taxa_order
    get_taxa_order <- function (tree_view) {
      df <- tree_view$data
      with(df, {
        i = order(y, decreasing = T)
        label[i][isTip[i]]
      })
    }
    
    # tips on tree plot from top
    tips_from_top <- get_taxa_order(p)
    
    # define tipdf in same order as tips on tree
    tipdf <- tree_tbl[idx_x, ]
    idx_top <- unname(sapply(tips_from_top, function(x) {
      which(tipdf$label == x)
    }))
    tipdf <- tipdf[idx_top, ]
    tipdf$label <- factor(tipdf$label, levels = rev(tipdf$label))
    
    for (i in seq_along(heatmap_var)) {
      
      if (length(heatmap_offset) == 1) {
        heatmap_offset_i <- heatmap_offset
      } else if (length(heatmap_offset) == length(heatmap_var)) {
        heatmap_offset_i <- heatmap_offset[i]
      } else {
        stop(paste0(
          "'heatmap_offset' must have length one or the same length as ",
          "'heatmap_var'"
        ))
      }
      
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
        if (length(unique(hmdf$group)) == 1) {
          hmdf_colors <- "grey"
          names(hmdf_colors) <- unique(hmdf$group)
        } else {
        hmdf_colors <- qualpalr::qualpal(
          length(unique(hmdf$group)), colorspace = "pretty")$hex
        }

        names(hmdf_colors) <- unique(hmdf$group)
        hmcolors <- "set"
      }

      var <- heatmap_var[i]
      
      hm_list[[i]] <- ggplot(tipdf, aes(x = "", y = label)) + 
        geom_tile(aes(fill = .data[[var]]))
      
      if (heatmap_text == TRUE) {
        hm_list[[i]] <- hm_list[[i]] +
          geom_text(aes(label = .data[[var]]), size = tiplab_size)
      }
      
      hm_list[[i]] <- hm_list[[i]] +
        xlab("") +
        ylab("") +
        theme(
          panel.background = element_rect(fill = "white"),
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
        )
      
      if (heatmap_var[i] %in% names(legend_breaks)) {
        breaks <- legend_breaks[[heatmap_var[i]]]
      } else {
        breaks <- names(hmdf_colors)
      }
      
      if (heatmap_var[i] %in% names(legend_guides)) {
        hm_list[[i]] <- hm_list[[i]] + scale_fill_manual(
          name = heatmap_var[i],
          values = hmdf_colors,
          breaks = breaks,
          na.translate = FALSE,
          guide = rlang::exec(guide_legend, !!!legend_guides[[heatmap_var[i]]])
        )
      } else {
        hm_list[[i]] <- hm_list[[i]] + scale_fill_manual(
          name = heatmap_var[i],
          values = hmdf_colors,
          breaks = breaks,
          na.translate = FALSE
        )
      }
      
    }
    
    p <- p + hm_list +
      patchwork::plot_layout(
        nrow = 1,
        widths = c(20, rep(heatmap_width, times = length(heatmap_var))),
        guides = "collect"
      )

  }
  p <- p &
    theme(
      plot.margin = unit(c(0, 0, 0, heatmap_offset), "points"),
      legend.position = legend_position,
      legend.box = legend_box
    )
  return(p)
}
