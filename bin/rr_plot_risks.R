library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-c", "--countlist"),
    type = "character",
    help = "Path to an rds file containing the counts for risk analysis."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    countlist = "countlist.rds"
  )
}

library(ggplot2)

countlist <- readRDS(args$countlist)

countdf <- countlist$countdf
mrca_cat_char = countlist$mrca_cat_char
geodist_threshold <- countlist$geodist_threshold

mrca_cat_char <- unique(countdf$mrca)

# convert counts to probabilities

probdf_type1 <- dplyr::bind_cols(
  countdf[,c(
    "subsample",
    "bootstrap",
    "mrca"
  )],
  as.data.frame(
    t(
      apply(countdf[,c(
        "same_country",
        "neighbors",
        "not_neighbors",
        "different_continent")], 1, function(x) {
          signif(x/sum(x), 4)
        })
    )
  )
)

probdf_type2 <- dplyr::bind_cols(
  countdf[,c(
    "subsample",
    "bootstrap",
    "mrca"
  )],
  as.data.frame(
    t(
      apply(countdf[,c(
        "same_country",
        "close_countries",
        "distant_countries",
        "different_continent")], 1, function(x) {
          signif(x/sum(x), 4)
        })
    )
  )
)

probdf_type3 <- dplyr::bind_cols(
  countdf[,c(
    "subsample",
    "bootstrap",
    "mrca"
  )],
  as.data.frame(
    t(
      apply(countdf[,c(
        "same_city",
        "different_city",
        "different_country")], 1, function(x) {
          signif(x/sum(x), 4)
        })
    )
  )
)

# convert probabilities to relative risks

rrdf_type1 <- dplyr::bind_cols(
  countdf[,c(
    "subsample",
    "bootstrap",
    "mrca"
  )],
  as.data.frame(
    t(
      apply(probdf_type1[,c(
        "same_country",
        "neighbors",
        "not_neighbors",
        "different_continent")], 1, function(x) {
          # this is where we hardcode "not neighbors" as the reference
          signif(x/x[3], 4)
        })
    )
  )
)

rrdf_type2 <- dplyr::bind_cols(
  countdf[,c(
    "subsample",
    "bootstrap",
    "mrca"
  )],
  as.data.frame(
    t(
      apply(probdf_type2[,c(
        "same_country",
        "close_countries",
        "distant_countries",
        "different_continent")], 1, function(x) {
          # this is where we hardcode "distant_countries" as the reference
          signif(x/x[3], 4)
        })
    )
  )
)

rrdf_type3 <- dplyr::bind_cols(
  countdf[,c(
    "subsample",
    "bootstrap",
    "mrca"
  )],
  as.data.frame(
    t(
      apply(probdf_type3[,c(
        "same_city",
        "different_city",
        "different_country")], 1, function(x) {
          # this is where we hardcode "different_country" as the reference
          signif(x/x[3], 4)
        })
    )
  )
)

rrdf_type1_long <- tidyr::pivot_longer(
  rrdf_type1,
  cols = same_country:different_continent,
  names_to = "geo",
  values_to = "rr"
)
rrdf_type1_long$geo <- factor(rrdf_type1_long$geo, levels = c(
  "same_country", "neighbors", "not_neighbors", "different_continent"
))
rrdf_type1_long$mrca <- factor(rrdf_type1_long$mrca, levels = mrca_cat_char)
rrdf_type1_long$rr <- ifelse(
  rrdf_type1_long$rr %in% c(NA, NaN, Inf, -Inf),
  NA,
  rrdf_type1_long$rr)


rrdf_type2_long <- tidyr::pivot_longer(
  rrdf_type2,
  cols = same_country:different_continent,
  names_to = "geo",
  values_to = "rr"
)
rrdf_type2_long$geo <- factor(rrdf_type2_long$geo, levels = c(
  "same_country", "close_countries", "distant_countries", "different_continent"
))
rrdf_type2_long$mrca <- factor(rrdf_type2_long$mrca, levels = mrca_cat_char)
rrdf_type2_long$rr <- ifelse(
  rrdf_type2_long$rr %in% c(NA, NaN, Inf, -Inf),
  NA,
  rrdf_type2_long$rr)

rrdf_type3_long <- tidyr::pivot_longer(
  rrdf_type3,
  cols = same_city:different_country,
  names_to = "geo",
  values_to = "rr"
)
rrdf_type3_long$geo <- factor(rrdf_type3_long$geo, levels = c(
  "same_city", "different_city", "different_country"
))
rrdf_type3_long$mrca <- factor(rrdf_type3_long$mrca, levels = mrca_cat_char)
rrdf_type3_long$rr <- ifelse(
  rrdf_type3_long$rr %in% c(NA, NaN, Inf, -Inf),
  NA,
  rrdf_type3_long$rr)

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

plot_rr <- function(df) {
  ggplot(df, aes(geo, rr)) + 
    stat_summary(
      geom = "point", 
      fun.data = point_and_whiskers,
      size = 0.3
    ) +
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

g1 <- plot_rr(rrdf_type1_long) +
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

close_countries_label <- paste0("Between countries \n <",geodist_threshold, "km")
distant_countries_label <- paste0("Between countries \n >",geodist_threshold, "km (ref)")

g2 <- plot_rr(rrdf_type2_long) +
  scale_x_discrete(labels = c(
    "same_country" = "Within \n countries",
    "close_countries" = close_countries_label,
    "distant_countries" = distant_countries_label,
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

g3 <- plot_rr(rrdf_type3_long) +
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

