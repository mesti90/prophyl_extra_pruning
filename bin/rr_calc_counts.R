library(optparse)
rm(list = ls())

args_list <- list(
  make_option(
    c("-p", "--project_dir"),
    type = "character",
    help = "Path to project directory."
  ),
  make_option(
    c("-A", "--assemblies"),
    type = "character",
    help = "Path to the assemblies file."
  ),
  make_option(
    c("-a", "--assemblies_collapsed"),
    type = "character",
    help = "Path to the assemblies file used for subsampling."
  ),
  make_option(
    c("-c", "--colldist"),
    type = "character",
    help = "Matrix of pairs, distances between sample collection dates."
  ),
  make_option(
    c("-y", "--same_city"),
    type = "character",
    help = "Matrix of pairs, are samples from the same city."
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
    c("-D", "--phylodist_all"),
    type = "character",
    help = "Matrices of pairs, phylogenetic distances for ALL tips."
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
  ),
  make_option(
    c("-m", "--mrca_categories"),
    type = "character",
    help = "MRCA categories."
  ),
  make_option(
    c("-C", "--colldist_max"),
    type = "character",
    help = "Maximum collection date distance between samples."
  )
)

args_parser  <- OptionParser(option_list = args_list)

if (!interactive()) {
  args  <- parse_args(args_parser)
} else {
  args <- list(
    project_dir = "aci/prophyl",
    assemblies = "assemblies.tsv",
    assemblies_collapsed = "assemblies_for_country_rr.tsv",
    colldist = "colldist.rds",
    same_city = "same_city.rds",
    same_country = "same_country.rds",
    neighbors = "neighbors.rds",
    same_continent = "same_continent.rds",
    geodist = "geodist.rds",
    phylodist_all = "phylodist.rds",
    phylodist = "phylodist_list.rds",
    nboot = 1,
    mrca_categories = "0,6,12,40",
    colldist_max = 2
  )
}

library(devtools)
load_all(args$project_dir)

assemblies <- read.csv(args$assemblies, sep = "\t")
colldist <- readRDS(args$colldist)
same_city <- readRDS(args$same_city)
same_country <- readRDS(args$same_country)
neighbors <- readRDS(args$neighbors)
same_continent <- readRDS(args$same_continent)
geodist <- readRDS(args$geodist)
phylodist_all <- readRDS(args$phylodist_all)
phylodist_list <- readRDS(args$phylodist)
nboot <- as.numeric(args$nboot)

# MRCA windows on which to compute the relative risk
# This will define categories on the risk plot
mrca_categories <- args$mrca_categories
mrca_categories <- as.numeric(strsplit(mrca_categories, ", *")[[1]])

mrca_cat_char <- paste0(
  "(",
  mrca_categories[-length(mrca_categories)],
  ",",
  mrca_categories[-1],
  "]"
)

# include lowest
mrca_cat_char[1] <- gsub("^\\(", "[", mrca_cat_char[1])

foo <- function(phd, categories) {
  out <- cut(phd, breaks = categories, include.lowest = TRUE)
  out <- matrix(out, ncol = ncol(phd))
  row.names(out) <- row.names(phd)
  colnames(out) <- colnames(phd)
  return(out)
}

# convert phylogenetic distances to MRCA categories
phylodist_list <- lapply(phylodist_list, function(x) {
  foo(x, categories = mrca_categories)
})

# Maximum collection date distance between samples
colldist_max <- as.numeric(args$colldist_max)
# Convert collection date distance matrix to boolean matrix
colldist <- 1 * (colldist <= colldist_max)

# Threshold between close and distant countries
geodist_threshold <- 1000
# Close countries are different countries within threshold
geodist_close <- 1 * (geodist < geodist_threshold)
close_countries <- (1-same_country) * geodist_close
# Distant countries are different countries beyond threshold
distant_countries <- (1-same_country) * (1-geodist_close)

countdf_all <- data.frame()

########## TODO DRY WRAP THIS INTO A FUNCTION?

# convert phylogenetic distances to MRCA categories

# filter all phylo distances to only those in the collapsed assemblies
assemblies_collapsed <- read.csv(args$assemblies_collapsed, sep = "\t")
index <- which(colnames(phylodist_all) %in% assemblies_collapsed$assembly)
phylodist_all <- phylodist_all[index, index]

phylodist_all <- foo(phylodist_all, mrca_categories)

# shrink comparison matrices to phylodist rows and columns
colldist_sub <- shrink_matrix(colldist, phylodist_all)
same_city_sub <- shrink_matrix(same_city, phylodist_all)
same_country_sub <- shrink_matrix(same_country, phylodist_all)
neighbors_sub <- shrink_matrix(neighbors, phylodist_all)
close_countries_sub <- shrink_matrix(close_countries, phylodist_all)
distant_countries_sub <- shrink_matrix(distant_countries, phylodist_all)
same_continent_sub <- shrink_matrix(same_continent, phylodist_all)

# for each MRCA category
for (k in mrca_cat_char) {
  phylodist_k <- 1 * (phylodist_all == k)
  ## CALCULATE BOOLEAN MATRICES
  # same city
  same_city_sub2 <- same_city_sub *
    colldist_sub * # within colldist timeframe
    phylodist_k # within MRCA range
  # different city, same country
  different_city_sub2 <- (1-same_city_sub) *
    same_country_sub *
    colldist_sub * # within colldist timeframe
    phylodist_k # within MRCA range
  # same country
  same_country_sub2 <- same_country_sub *
    colldist_sub * # within colldist timeframe
    phylodist_k # within MRCA range
  # different country, same continent
  different_country_same_continent_sub2 <- (1-same_country_sub2) *
    same_continent_sub *
    colldist_sub * # within colldist timeframe
    phylodist_k # within MRCA range
  # different country, irrespective of continent
  different_country_sub2 <- (1-same_country_sub2) *
    colldist_sub * # within colldist timeframe
    phylodist_k # within MRCA range
  # neighbors, same continent
  neighbors_sub2 <- (1-same_country_sub) *
    neighbors_sub *
    same_continent_sub *
    colldist_sub * # within colldist timeframe
    phylodist_k # within MRCA range
  # different countries, not neighbors, same continent
  not_neighbors_sub2 <- (1-same_country_sub) *
    (1-neighbors_sub) *
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
  different_continent_sub2 <- (1-same_continent_sub) * 
    colldist_sub * # within colldist timeframe
    phylodist_k # within MRCA range

  ddf <- data.frame(
    same_city = sum(same_city_sub2, na.rm = TRUE) / 2,
    different_city = sum(different_city_sub2, na.rm = TRUE) / 2,
    same_country = sum(same_country_sub2, na.rm = TRUE) / 2,
    different_country_same_continent = sum(different_country_same_continent_sub2, na.rm = TRUE) / 2,
    different_country = sum(different_country_sub2, na.rm = TRUE) / 2,
    neighbors = sum(neighbors_sub2, na.rm = TRUE) / 2,
    not_neighbors = sum(not_neighbors_sub2, na.rm = TRUE) / 2,
    close_countries = sum(close_countries_sub2, na.rm = TRUE) / 2,
    distant_countries = sum(distant_countries_sub2, na.rm = TRUE) / 2,
    different_continent = sum(different_continent_sub2, na.rm = TRUE) / 2
  )

  newdf <- dplyr::bind_cols(
    data.frame(
      subsample = "All",
      bootstrap = NA,
      mrca = k
    ),
    ddf
  )
      
  countdf_all <- dplyr::bind_rows(
    countdf_all,
    newdf
  )
}

countlist_all <- list(
  countdf = countdf_all,
  mrca_cat_char = mrca_cat_char,
  geodist_threshold = geodist_threshold
)

saveRDS(countlist_all, file = "countlist_all.rds")

##################################

countdf <- data.frame()
# for each subsample
for (i in seq_along(phylodist_list)) {
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
    same_city_sub <- shrink_matrix(same_city, phylodist)
    same_country_sub <- shrink_matrix(same_country, phylodist)
    neighbors_sub <- shrink_matrix(neighbors, phylodist)
    close_countries_sub <- shrink_matrix(close_countries, phylodist)
    distant_countries_sub <- shrink_matrix(distant_countries, phylodist)
    same_continent_sub <- shrink_matrix(same_continent, phylodist)
    
    # for each MRCA category
    for (k in mrca_cat_char) {
      phylodist_k <- 1 * (phylodist == k)
      ## CALCULATE BOOLEAN MATRICES
      # same city
      same_city_sub2 <- same_city_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # different city, same country
      different_city_sub2 <- (1-same_city_sub) *
        same_country_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # same country
      same_country_sub2 <- same_country_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # different country, same continent
        different_country_same_continent_sub2 <- (1-same_country_sub2) *
        same_continent_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # different country, irrespective of continent
      different_country_sub2 <- (1-same_country_sub2) *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # neighbors, same continent
      neighbors_sub2 <- (1-same_country_sub) *
        neighbors_sub *
        same_continent_sub *
        colldist_sub * # within colldist timeframe
        phylodist_k # within MRCA range
      # different countries, not neighbors, same continent
      not_neighbors_sub2 <- (1-same_country_sub) *
        (1-neighbors_sub) *
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
      different_continent_sub2 <- (1-same_continent_sub) * 
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
      
      subset_matrix <- function(matrix, value) {
        index_x <- vector()
        index_y <- vector()
        for (i in 1:nrow(matrix)) {
          for (j in 1:ncol(matrix)) {
            if (!is.na(matrix[i,j]) & matrix[i,j] == value) {
              index_x <- c(index_x, i)
              index_y <- c(index_y, j)
            }
          }
        }
        subset_mat <- matrix[index_x, index_y]
        return(subset_mat)
      }
      
      # type 1 - same country, neighbor, not neighbor, different continent
      # testthat::expect_false(matrix_overlap(
      #   same_country_sub2,
      #   neighbors_sub2
      # ))
      # testthat::expect_false(matrix_overlap(
      #   same_country_sub2,
      #   not_neighbors_sub2
      # ))
      # testthat::expect_false(matrix_overlap(
      #   same_country_sub2,
      #   different_continent_sub2
      # ))
      # testthat::expect_false(matrix_overlap(
      #   neighbors_sub2,
      #   not_neighbors_sub2
      # ))
      # testthat::expect_false(matrix_overlap(
      #   neighbors_sub2,
      #   different_continent_sub2
      # ))
      # testthat::expect_false(matrix_overlap(
      #   not_neighbors_sub2,
      #   different_continent_sub2
      # ))
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

      # type3 - same city, different city, different country
      testthat::expect_false(matrix_overlap(
        same_city_sub2,
        different_city_sub2,
        different_country_sub2
      ))

      ddf <- data.frame(
        same_city = sum(same_city_sub2, na.rm = TRUE) / 2,
        different_city = sum(different_city_sub2, na.rm = TRUE) / 2,
        same_country = sum(same_country_sub2, na.rm = TRUE) / 2,
        different_country_same_continent = sum(different_country_same_continent_sub2, na.rm = TRUE) / 2,
        different_country = sum(different_country_sub2, na.rm = TRUE) / 2,
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
        ddf$different_country_same_continent,
        ddf$different_continent
      )

      type3_sumcount <- sum(
        ddf$same_city,
        ddf$different_city,
        ddf$different_country
      )
      
      testthat::expect_true(type1_sumcount == type2_sumcount)
      
      newdf <- dplyr::bind_cols(
        data.frame(
          subsample = names(phylodist_list[i]),
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

countdf <- countdf[order(countdf$subsample),]

countlist <- list(
  countdf = countdf,
  mrca_cat_char = mrca_cat_char,
  geodist_threshold = geodist_threshold
)

saveRDS(countlist, file = "countlist.rds")
