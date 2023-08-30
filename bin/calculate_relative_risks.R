library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("-a", "--assemblies"),
    type = "character",
    help = "Path to the assemblies file."
  ),
  make_option(
    c("-c", "--colldist"),
    type = "character",
    help = "Matrix of pairs, distances between sample collection dates."
  ),
  make_option(
    c("-r", "--same_country"),
    type = "character",
    help = "Matrix of pairs, are samples from the same country."
  ),
  make_option(
    c("-n", "--neighbors"),
    type = "character",
    help = "Matrix of pairs, are samples from the neighboring countries."
  ),
  make_option(
    c("-t", "--same_continent"),
    type = "character",
    help = "Matrix of pairs, are samples from the same continent."
  ),
  make_option(
    c("-g", "--geodist"),
    type = "character",
    help = "Matrix of pairs, geographical distances between samples."
  ),
  make_option(
    c("-d", "--phylodist"),
    type = "character",
    help = "List of matrices, matrices of pairs, phylogenetic distances."
  ),
  make_option(
    c("-b", "--nboot"),
    type = "character",
    help = "Number of bootstrap replicates for each simulated tree."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "~/Methods/prophyl",
    assemblies = "assemblies.tsv",
    colldist = "colldist.rds",
    same_country = "same_country.rds",
    neighbors = "neighbors.rds",
    same_continent = "same_continent.rds",
    geodist = "geodist.rds",
    phylodist = "phylodist_list.rds",
    nboot = 1
  )
}

library(devtools)
load_all(args$project_dir)

assemblies <- read.csv(args$assemblies, sep = "\t")
colldist <- readRDS(args$colldist)
same_country <- readRDS(args$same_country)
neighbors <- readRDS(args$neighbors)
same_continent <- readRDS(args$same_continent)
geodist <- readRDS(args$geodist)
phylodist_list <- readRDS(args$phylodist)
nboot <- as.numeric(args$nboot)

# MRCA windows on which to compute the relative risk
# This will define categories on the risk plot
mrca_categories <- c(0, 4, 8, 16, Inf)

mrca_cat_char <- paste0(
  "(",
  mrca_categories[-length(mrca_categories)],
  ",",
  mrca_categories[-1],
  "]")

# convert phylogenetic distances to MRCA categories
phylodist_list <- lapply(phylodist_list, function(x) {
  out <- cut(x, breaks = mrca_categories)
  out <- matrix(out, ncol = ncol(x))
  row.names(out) <- row.names(x)
  colnames(out) <- colnames(x)
  return(out)
})

# Maximum collection date distance between samples
colldist_max <- 2
# Convert collection date distance matrix to boolean matrix
colldist <- 1 * (colldist < colldist_max)

# Threshold between close and distant countries
geodist_threshold <- 1000
# Close countries are different countries within threshold
geodist_close <- 1 * (geodist < geodist_threshold)
close_countries <- (1*!same_country) * geodist_close
# Distant countries are different countries beyond threshold
distant_countries <- (1*!same_country) * (1*!geodist_close)

# Define function for shrinking a matrix to have the same rows and columns as
# the smaller matrix
shrink_matrix <- function(matrix, small_matrix) {
  index <- unname(sapply(colnames(small_matrix), function(x) {
    which(colnames(matrix) == x)
  }))
  shrinked_matrix <- matrix[index, index]
  return(shrinked_matrix)
}

countdf <- data.frame()
# for each subsample
for (i in 1:length(phylodist_list)) {
  # for each bootstrap replicate
  for (j in 1:nboot) {
    # use the respective phylodist matrix or generate a bootstrapped matrix
    if (nboot == 1) {
      phylodist <- phylodist_list[[i]]
    } else {
      index <- sample(
        1:ncol(phylodist_list[[i]]),
        ncol(phylodist_list[[i]]),
        replace = TRUE
      )
      phylodist <- phylodist_list[[i]][index,index]
    }
    # shrink comparison matrices to phylodist rows and columns
    colldist_sub <- shrink_matrix(colldist, phylodist)
    same_country_sub <- shrink_matrix(same_country, phylodist)
    neighbors_sub <- shrink_matrix(neighbors, phylodist)
    close_countries_sub <- shrink_matrix(close_countries, phylodist)
    distant_countries_sub <- shrink_matrix(distant_countries, phylodist)
    same_continent_sub <- shrink_matrix(same_continent, phylodist)
    
    # for each MRCA category
    for (k in mrca_cat_char) {
      phylodist_k <- 1 * (phylodist == k)
      ## CALCULATE BOOLEAN MATRICES
      # same country
      same_country_sub2 <- same_country_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # neighbors, same continent
      neighbors_sub2 <- (1*!same_country_sub) *
        neighbors_sub *
        same_continent_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # different countries, not neighbors, same continent
      not_neighbors_sub2 <- (1*!same_country_sub) *
        (1*!neighbors_sub2) *
        same_continent_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # close_countries, same continent
      close_countries_sub2 <- close_countries_sub *
        same_continent_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # distant countries, same continent
      distant_countries_sub2 <- distant_countries_sub * 
        same_continent_sub *
        colldist_sub *  # within colldist timeframe
        phylodist_k # within MRCA range
      # different continent
      different_continent_sub2 <- (1*!same_continent_sub) * 
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      
      # consistency checks - no overlaps between exlusive categories
      
      matrix_overlap <- function(...) {
        args <- list(...)
        if (length(args) < 2) {
          stop("At least two arguments are required")
        }
        smat <- args[[1]]
        for (i in 2:length(args)) {
          smat <- smat + args[[i]]
        }
        smat_tab <- table(smat)
        all(names(smat_tab) %in% c("0", "1")) == FALSE
      }
      
      # type 1 - same country, neighbor, not neighbor, different continent
      testthat::expect_false(matrix_overlap(
        same_country_sub2,
        neighbors_sub2
      ))
      testthat::expect_false(matrix_overlap(
        same_country_sub2,
        not_neighbors_sub2
      ))
      testthat::expect_false(matrix_overlap(
        same_country_sub2,
        different_continent_sub2
      ))
      testthat::expect_false(matrix_overlap(
        neighbors_sub2,
        not_neighbors_sub2
      ))
      testthat::expect_false(matrix_overlap(
        neighbors_sub2,
        different_continent_sub2
      ))
      testthat::expect_false(matrix_overlap(
        not_neighbors_sub2,
        different_continent_sub2
      ))
      testthat::expect_false(matrix_overlap(
        same_country_sub2,
        neighbors_sub2,
        not_neighbors_sub2,
        different_continent_sub2
      ))
      # type 2- same country, close country, distant county, different continent
      testthat::expect_false(matrix_overlap(
        same_country_sub2,
        close_countries_sub2,
        distant_countries_sub2,
        different_continent_sub2
      ))
      
      ddf <- data.frame(
        same_country = sum(same_country_sub2, na.rm = TRUE) / 2,
        neighbors = sum(neighbors_sub2, na.rm = TRUE) / 2,
        not_neighbors = sum(not_neighbors_sub2, na.rm = TRUE) / 2,
        close_countries = sum(close_countries_sub2, na.rm = TRUE) / 2,
        distant_countries = sum(distant_countries_sub2, na.rm = TRUE) / 2,
        different_continent = sum(different_continent_sub2, na.rm = TRUE) / 2
      )
      
      # consistency check - each type covers all pairs
      type1_sumcount <- sum(
        ddf$same_country, 
        ddf$neighbors,
        ddf$not_neighbors,
        ddf$different_continent
      )
      
      type2_sumcount <- sum(
        ddf$same_country,
        ddf$close_countries,
        ddf$distant_countries,
        ddf$different_continent
      )
      
      testthat::expect_true(type1_sumcount == type2_sumcount)
      
      # fill up rows of the data frame
      # TODO subsample should contain the real name of the subsample
      # it is unclear in what order phylodist list contains the subsamples.
      newdf <- dplyr::bind_cols(
        data.frame(
          subsample = i,
          bootstrap = j,
          mrca = k
        ),
        ddf
      )
      
      countdf <- dplyr::bind_rows(
        countdf,
        newdf
      )
    }
  }
}

saveRDS(countdf, file = "counts.rds")

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
