library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-c", "--countlist"),
    type = "character",
    help = "Path to an rds file containing the counts for risk analysis."
  ),
  make_option(
    c("-C", "--countlist_all"),
    type = "character",
    help = "Path to an rds file containing the counts for risk analysis."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    countlist = "countlist.rds",
    countlist_all = "countlist_all.rds"
  )
}

library(ggplot2)

countlist <- readRDS(args$countlist)
countlist_all <- readRDS(args$countlist_all)

countdf <- countlist$countdf
mrca_cat_char = countlist$mrca_cat_char
geodist_threshold <- countlist$geodist_threshold

mrca_cat_char <- unique(countdf$mrca)

# convert counts to probabilities

get_probdf <- function(df, vars) {
  out <- dplyr::bind_cols(
    df[,c(
      "subsample",
      "bootstrap",
      "mrca"
    )],
    as.data.frame(
      t(
        apply(df[,vars], 1, function(x) {
          signif(x/sum(x), 4)
        })
      )
    )
  )
  return(out)
}

type1_vars <- c(
  "same_country",
  "neighbors",
  "not_neighbors",
  "different_continent"
)

type2_vars <- c(
  "same_country",
  "different_country_same_continent",
  "different_continent"
)

type3_vars <- c(
  "same_city",
  "different_city",
  "different_country"
)

probdf_type1 <- get_probdf(countdf, type1_vars)
probdf_type1_all <- get_probdf(countlist_all$countdf, type1_vars)

probdf_type2 <- get_probdf(countdf, type2_vars)
probdf_type2_all <- get_probdf(countlist_all$countdf, type2_vars)

probdf_type3 <- get_probdf(countdf, type3_vars)
probdf_type3_all <- get_probdf(countlist_all$countdf, type3_vars)

# convert probabilities to relative risks

get_rrdf <- function(df, vars, ref) {
  index_ref <- which(names(df[, vars]) == ref)
  out <- dplyr::bind_cols(
    df[,c(
      "subsample",
      "bootstrap",
      "mrca"
    )],
    as.data.frame(
      t(
        apply(df[,vars], 1, function(x) {
          signif(x/x[index_ref], 4)
        })
      )
    )
  )
  return(out)
}

rrdf_type1 <- get_rrdf(probdf_type1, type1_vars, "not_neighbors")
rrdf_type1_all <- get_rrdf(probdf_type1_all, type1_vars, "not_neighbors")

rrdf_type2 <- get_rrdf(probdf_type2, type2_vars, "different_country_same_continent")
rrdf_type2_all <- get_rrdf(probdf_type2_all, type2_vars, "different_country_same_continent")

rrdf_type3 <- get_rrdf(probdf_type3, type3_vars, "different_country")
rrdf_type3_all <- get_rrdf(probdf_type3_all, type3_vars, "different_country")

format_long <- function(df, vars) {
  cols <- names(df)[which(names(df) %in% vars)]
  long <- tidyr::pivot_longer(
    df,
    cols = cols,
    names_to = "geo",
    values_to = "rr"
  )
  long$geo <- factor(long$geo, levels = vars)
  long$mrca <- factor(long$mrca, levels = mrca_cat_char)
  long$rr <- ifelse(
    long$rr %in% c(NA, NaN, Inf, -Inf),
    NA,
    long$rr
  )
  return(long)
}

rrdf_type1_long <- format_long(rrdf_type1, type1_vars)
rrdf_type1_long_all <- format_long(rrdf_type1_all, type1_vars)

rrdf_type2_long <- format_long(rrdf_type2, type2_vars)
rrdf_type2_long_all <- format_long(rrdf_type2_all, type2_vars)

rrdf_type3_long <- format_long(rrdf_type3, type3_vars)
rrdf_type3_long_all <- format_long(rrdf_type3_all, type3_vars)

point_and_whiskers <- function(x) {
  y <- median(x)
  ymin <- quantile(x, 0.05)
  ymax <- quantile(x, 0.95)
  return(data.frame(
    "y" = y,
    "ymin" = ymin,
    "ymax" = ymax
  ))
}

plot_rr <- function(df, df_all) {
  ggplot(df, aes(geo, rr)) + 
    stat_summary(
      geom = "point", 
      fun.data = point_and_whiskers,
      size = 0.3
    ) + 
    # uncomment this to include risk values for the full data set
    # stat_summary(
    #   geom = "point", 
    #   fun.data = point_and_whiskers,
    #   size = 0.3, 
    #   data = df_all, 
    #   col = "#FF0000"
    # ) +
    stat_summary(
      geom = "errorbar",
      fun.data = point_and_whiskers,
      width = 0.1,
      linewidth = 0.1,
      col = "#000000"
    ) +
    geom_hline(yintercept = 1, col = "#FF0000", linewidth = 0.1) +
    facet_grid(mrca~.) +
    scale_y_log10() +
    ylab("Relative Risk") +
    xlab("") +
    theme_minimal() +
    theme(
      panel.background = element_rect(
        colour = "#000000",
        linewidth = 0.1
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(family = "helvetica", size = 5),
      axis.text = element_text(family = "helvetica", size = 5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.ticks = element_line(colour = "#000000", linewidth = 0.1),
      axis.ticks.length = unit(0.05, "cm"),
      strip.text = element_text(family = "helvetica", size = 5)
    ) 
}

g1 <- plot_rr(rrdf_type1_long, rrdf_type1_long_all) +
  scale_x_discrete(labels = c(
    "same_country" = "Within \n countries",
    "neighbors" = "Between \n neighbors",
    "not_neighbors" = "Between \n non-neighbors (ref)",
    "different_continent" = "Between \n continents"
  ))

try_g1 <- try(print(g1), silent = TRUE)
if (inherits(try_g1, "try-error")) {
  g1 <- ggplot()
}

ggsave(
  filename = "relative_risks_type1.pdf",
  plot = g1,
  units = "cm",
  width = 8,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = "relative_risks_type1.png",
  plot = g1,
  units = "cm",
  width = 8,
  height = 6
)

saveRDS(g1, "relative_risks_type1.rds")

between_countries_label <- paste0("Between countries")

g2 <- plot_rr(rrdf_type2_long, rrdf_type2_long_all) +
  scale_x_discrete(labels = c(
    "same_country" = "Within \n countries",
    "different_country_same_continent" = "Between \n countries",
    "different_continent" = "Between \n continents"
  ))

try_g2 <- try(print(g2), silent = TRUE)
if (inherits(try_g2, "try-error")) {
  g2 <- ggplot()
}

ggsave(
  filename = "relative_risks_type2.pdf",
  plot = g2,
  units = "cm",
  width = 8,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = "relative_risks_type2.png",
  plot = g2,
  units = "cm",
  width = 8,
  height = 6
)

saveRDS(g2, "relative_risks_type2.rds")

g3 <- plot_rr(rrdf_type3_long, rrdf_type3_long_all) +
  scale_x_discrete(labels = c(
    "same_country" = "Within \n countries",
    "close_countries" = close_countries_label,
    "distant_countries" = distant_countries_label,
    "different_continent" = "Between \n continents"
  ))

try_g3 <- try(print(g3), silent = TRUE)
if (inherits(try_g3, "try-error")) {
  g3 <- ggplot()
}

ggsave(
  filename = "relative_risks_type3.pdf",
  plot = g3,
  units = "cm",
  width = 8,
  height = 6,
  device = cairo_pdf
)

ggsave(
  filename = "relative_risks_type3.png",
  plot = g3,
  units = "cm",
  width = 8,
  height = 6
)

saveRDS(g3, "relative_risks_type3.rds")

