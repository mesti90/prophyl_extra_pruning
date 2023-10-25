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

library(dplyr)
library(ggplot2)

prettify_mrca <- function(mrca) {
  foo <- function(x) {
    sp <- strsplit(x, "\\[|\\(|,|\\)|\\]")[[1]]
    paste0(sp[2], " < MRCA (years) < ", sp[3])
  }
  unname(sapply(mrca, foo))
}

# import countlist
countlist <- readRDS(args$countlist)
# make mrca categories pretty
countlist$mrca_cat_char <- prettify_mrca(countlist$mrca_cat_char)
countlist$countdf$mrca <- prettify_mrca(countlist$countdf$mrca)

# import countlist for all tips
countlist_all <- readRDS(args$countlist_all)
countlist_all$mrca_cat_char <- prettify_mrca(countlist_all$mrca_cat_char)
countlist_all$countdf$mrca <- prettify_mrca(countlist_all$countdf$mrca)

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
  "different_city"
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

rrdf_type3 <- get_rrdf(probdf_type3, type3_vars, "different_city")
rrdf_type3_all <- get_rrdf(probdf_type3_all, type3_vars, "different_city")

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
  if (all(is.na(x))) {
    df <- data.frame(
      "y" = NA,
      "ymin" = NA,
      "ymax" = NA
    )
  } else {
    df <- data.frame(
      "y" = median(x, na.rm = TRUE),
      "ymin" = quantile(x, 0.05, na.rm = TRUE),
      "ymax" = quantile(x, 0.95, na.rm = TRUE)
    )
  }
  return(df)
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
    facet_wrap(mrca~., nrow = length(unique(df$mrca))) +
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
      axis.title = element_text(size = 5),
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.ticks = element_line(colour = "#000000", linewidth = 0.1),
      axis.ticks.length = unit(0.05, "cm"),
      strip.text = element_text(
        size = 5,
        margin = margin(0,0,0,0, "cm")
      )
    ) 
}

get_p_values <- function(gg, type) {
  out <- gg$data %>% 
    group_by(mrca, geo) %>% 
    summarise(
      higher = mean(rr > 1, na.rm = TRUE),
      lower = mean(rr < 1, na.rm = TRUE),
      .groups = "drop"
    )
  out <- dplyr::bind_cols(type = type, out)
  return(out)
}

g1 <- plot_rr(rrdf_type1_long, rrdf_type1_long_all) +
  scale_x_discrete(labels = c(
    "same_country" = "Within \n countries",
    "neighbors" = "Between \n neighbors",
    "not_neighbors" = "Between \n non-neighbors (ref)",
    "different_continent" = "Between \n continents"
  ))

if (inherits(try(ggplot_build(g1), silent = TRUE), "try-error")) {
  g1 <- ggplot()
}

if (inherits(try(ggplot_build(g1), silent = TRUE), "try-error")) {
  ps <- data.frame()
} else {
  ps <- get_p_values(g1, "type1")
}

ggsave(
  filename = "relative_risks_type1.pdf",
  plot = g1,
  units = "cm",
  width = 8,
  height = 6
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

if (inherits(try(ggplot_build(g2), silent = TRUE), "try-error")) {
  g2 <- ggplot()
}

if (inherits(try(ggplot_build(g2), silent = TRUE), "try-error")) {
  ps <- dplyr::bind_rows(ps, data.frame())
} else {
  ps <- dplyr::bind_rows(ps, get_p_values(g2, "type2"))
}

ggsave(
  filename = "relative_risks_type2.pdf",
  plot = g2,
  units = "cm",
  width = 8,
  height = 6
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
    "same_city" = "Within \n cities",
    "different_city" = "Within \n countries"
  ))

if (inherits(try(ggplot_build(g3), silent = TRUE), "try-error")) {
  g3 <- ggplot()
}

if (inherits(try(ggplot_build(g3), silent = TRUE), "try-error")) {
  ps <- dplyr::bind_rows(ps, data.frame())
} else {
  ps <- dplyr::bind_rows(ps, get_p_values(g3, "type3"))
}

ggsave(
  filename = "relative_risks_type3.pdf",
  plot = g3,
  units = "cm",
  width = 8,
  height = 6
)

ggsave(
  filename = "relative_risks_type3.png",
  plot = g3,
  units = "cm",
  width = 8,
  height = 6
)

saveRDS(g3, "relative_risks_type3.rds")

write.table(
  ps, 
  file = "p_values.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)