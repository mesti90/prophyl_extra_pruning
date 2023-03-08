rm(list = ls())

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  df_path <- args[1]
  collapse_by <- args[2]
} else {
  df_path <- "~/Methods/prophyl-test/results/validate_input/assemblies.tsv"
  collapse_by <- "geo_date"
}

df <- read.csv(df_path, sep = "\t")

# add new variable for combined serotype
# TODO generalise to cases where there is no such variable
df$serotype <- paste(df$mlst, df$k_serotype, sep = "_")

set.seed(0)
if (collapse_by == "geo_date") {
  # add new temporal variables - yearweek
  df$yearweek <- paste(
    df$collection_year, lubridate::week(df$collection_day), sep = "_")
  df$yearweek <- ifelse(grepl("NA", df$yearweek), NA, df$yearweek)
  # add new temporal variables - yearmonth
  df$yearmonth <- paste(
    df$collection_year, lubridate::month(df$collection_day), sep = "_")
  df$yearmonth <- ifelse(grepl("NA", df$yearmonth), NA, df$yearmonth)
  # collapse samples where city is known
    cities <- unique(df$city)[which(!is.na(unique(df$city)))]
  keep_city <- vector()
  for (i in cities) {
    geo <- df[which(df$city == i),]
    serotypes <- unique(geo$serotype)
    for (j in serotypes) {
      geo_sero <- geo[which(geo$serotype == j),]
      for (k in unique(geo_sero$yearmonth)) {
        if (is.na(k)) {
          # if either year or month or both are not known, keep assembly
          keep_city <- c(
            keep_city, geo_sero$assembly[which(is.na(geo_sero$yearmonth))])
        } else {
          # otherwise keep one assembly per year per month
          geo_sero_yearmonth <- geo_sero[which(geo_sero$yearmonth == k),]
          keep_city <- c(
            keep_city, sample(geo_sero_yearmonth$assembly, 1))
        }
      }
    }
  }
  # collapse samples where city is not known but country is known
  countries <- unique(df$country)[which(!is.na(unique(df$country)))]
  keep_country <- vector()
  for (i in countries) {
    # only consider samples where the country is known but the city is not
    geo <- df[which(df$country == i & is.na(df$city)),]
    serotypes <- unique(geo$serotype)
    for (j in serotypes) {
      geo_sero <- geo[which(geo$serotype == j),]
      for (k in unique(geo_sero$yearweek)) {
        if (is.na(k)) {
          # if either year or week or both are not known, keep assembly
          keep_country <- c(
            keep_country, geo_sero$assembly[which(is.na(geo_sero$yearweek))])
        } else {
          # otherwise keep one assembly per year per week
          geo_sero_yearweek <- geo_sero[which(geo_sero$yearweek == k),]
          keep_country <- c(
            keep_country, sample(geo_sero_yearweek$assembly, 1))
        }
      }
    }
  }
  index <- which(
    df$assembly %in% keep_city |
    df$assembly %in% keep_country |
    is.na(df$country)
  )
  df <- df[index,]
}

if (!interactive()) {
  saveRDS(df, file = "assemblies_collapsed_outbreaks.rds")
}
