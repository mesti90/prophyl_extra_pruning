# Original script from:
# https://github.com/noemielefrancq/Global_spread_Listeria_monocytogenes_CC1

#######################################################
## Figure 4A: Relative risk by interval, across different location
#######################################################
### Author: Noemie Lefrancq
### Date creation: 03/02/2020
### Last modification: 17/10/2021
#######################################################

# This script has been modified to fit in the analysis pipeline

rm(list = ls())

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  data_seq_path <- args[1]
  time_mat_path <- args[2]
  geo_mat_country_path <- args[3]
  geo_mat_continent_path <- args[4]
  geo_mat_km_centroids_path <- args[5]
  sim_mats_path <- args[6]
  nboot <- as.numeric(args[7])
  project_dir <- args[8]
} else {
  test_dir <- "~/Methods/prophyl-tests/test-calculate_relative_risks"
  data_seq_path <- paste0(test_dir, "/assemblies.tsv")
  time_mat_path <- paste0(test_dir, "/results/calculate_distances/colldist.rds")
  geo_mat_country_path <- paste0(
    test_dir, "/results/calculate_distances/same_country.rds")
  geo_mat_continent_path <- paste0(
    test_dir, "/results/calculate_distances/same_continent.rds")
  geo_mat_km_centroids_path <- paste0(
    test_dir, "/results/calculate_distances/geodist.rds")
  sim_mats_path <- paste0(
    test_dir, "/results/calculate_distances/phylodist_list.rds")
  nboot <- 1
  project_dir <- "~/Methods/prophyl"
}

library(devtools)
load_all(project_dir)

#####################################################################
## Load relevant datasets and matrices
#####################################################################
## Metadata
data.seq = read.csv(data_seq_path, sep = "\t")
## Time
time_mat = readRDS(time_mat_path)
## Geography
geo_mat_country = readRDS(geo_mat_country_path)
geo_mat_continent = readRDS(geo_mat_continent_path)
geo_mat_km_centroids = readRDS(geo_mat_km_centroids_path)
## Genetic distances
sim.mats <- readRDS(sim_mats_path)
nsim = length(sim.mats)

#####################################################################
## Parameters for the computation of the relative risks
#####################################################################
## Number of bootstrap event to perform, of each tree
nboot = nboot

## MRCA windows on which to compute the relative risk
Pmax <- c(5, 10, 20, 40, 1000) ## max windows
Pmin <- c(0, 5, 10, 20, 40) ## min windows
pmid <- (Pmin+Pmax)/2 ## mid-point
int <- c(0, 5, 10, 20, 40, 1000)
l = length(Pmax) ## number of intervals

n_steps = 4 ## Number of location matrix to consider, here:
## 1- same country in EU
## 2- different countries <1000km in the same continent
## 3- different countries >1000 km in the same continent (reference)
## 4- different continents

## Set the boot matrix to save the results
rr = matrix(NA,l*n_steps, nboot*nsim)

#####################################################################
## Compute relative risks, for each location
#####################################################################
for (ii in 1:nsim) {
  # Choose isolates
  a <- unname(sapply(colnames(sim.mats[[ii]]), function(x) {
    which(data.seq$assembly == x)
  }))
  ## Time between isolates: max 2 years
  time_mat2 = time_mat[a,a]<=2  
  # Choose MRCA matrix
  MRCA_mat = sim.mats[[ii]]
  MRCA_mat2 <- MRCA_mat
  nseq = length(a)
  # Reference: different countries >1000 km in the same continent
  ref <- unname(sapply(colnames(sim.mats[[ii]]), function(x) {
    which(data.seq$assembly == x)
  }))
  geo_mat_ref = (1-geo_mat_country[ref,ref])*(geo_mat_km_centroids[a,a]>1000)*(geo_mat_continent[ref,ref])
  geo_mat_ref[which(geo_mat_ref == 0)] = NA
  ## Time between isolates: max 2 years
  time_mat2_ref = time_mat[ref,ref]<=2
  ## Choose MRCA matrix
  MRCA_mat2_ref = MRCA_mat
  nseq_ref = length(ref)
  
  for (j in 1:n_steps) {
    if (j == 1) {
      # same country
      geo_mat = geo_mat_country[a,a]
      ##Bootstrap to create the ci
      for (i in (1:nboot)){
        if (nboot == 1) {
          tmp <- 1:nseq
          tmp_ref <- 1:nseq_ref
        } else {
          tmp = sample(nseq, nseq, replace = T)
          tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
        }
        rr.out = ratio_bootstrap_dist_discrete_auto(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref, Pmax, Pmin)
        rr[((j*l-(l-1)):(j*l)),((ii-1)*nboot + i)] = rr.out
      }
    }
    if (j == 2) {
      # different countries in the same continent, less than 1000km apart
      geo_mat = (1-geo_mat_country[a,a])*(geo_mat_km_centroids[a,a]<=1000)*(geo_mat_continent[a,a])
      geo_mat[which(geo_mat == 0)] = NA
      ##Bootstrap to create the ci
      for (i in (1:nboot)){
        if (nboot == 1) {
          tmp <- 1:nseq
          tmp_ref <- 1:nseq_ref
        } else {
          tmp = sample(nseq, nseq, replace = T)
          tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
        }
        rr.out = ratio_bootstrap_dist_discrete_auto(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref, Pmax, Pmin)
        rr[((j*l-(l-1)):(j*l)),((ii-1)*nboot + i)] = rr.out
      }
    }
    if (j == 3) {
      # different countries in the same continent, more than 1000km apart (reference)
      geo_mat = (1-geo_mat_country[a,a])*(geo_mat_km_centroids[a,a]>1000)*(geo_mat_continent[a,a])
      geo_mat[which(geo_mat == 0)] = NA
      ##Bootstrap to create the ci
      for (i in (1:nboot)){
        if (nboot == 1) {
          tmp <- 1:nseq
          tmp_ref <- 1:nseq_ref
        } else {
          tmp = sample(nseq, nseq, replace = T)
          tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
        }
        rr.out = ratio_bootstrap_dist_discrete_auto(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref, Pmax, Pmin)
        rr[((j*l-(l-1)):(j*l)),((ii-1)*nboot + i)] = rr.out
      }
    }
    if (j == 4) {
      # different continents
      geo_mat = 1-geo_mat_continent[a,a]
      ##Bootstrap to create the ci
      for (i in (1:nboot)){
        if (nboot == 1) {
          tmp <- 1:nseq
          tmp_ref <- 1:nseq_ref
        } else {
          tmp = sample(nseq, nseq, replace = T)
          tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
        }
        rr.out = ratio_bootstrap_dist_discrete_auto(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref, Pmax, Pmin)
        rr[((j*l-(l-1)):(j*l)),((ii-1)*nboot + i)] = rr.out
      }
    }
    
  }
}

#####################################################################

#####################################################################
## Write results
#####################################################################
res = list('rr' = rr,
           'int' = int,
           'nboot' = nboot,
           'nsim' = nsim)
class(res) <- "rrlist"

if (!interactive()) {
  # export results
  saveRDS(res, "relative_risks.rds")
  # plot results - pdf
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  pdf(file = "relative_risks.pdf")
  plot_rr(res)
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  # plot results - png
  png(file = "relative_risks.png")
  plot_rr(res)
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
}
