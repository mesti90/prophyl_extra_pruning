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

args <- commandArgs(trailingOnly = TRUE)

## Main function
ratio_bootstrap_dist_discrete_auto <- function(x, y, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref){
  ## This function computes the relative risk, for a given 
  ##     -location matrix (geo_mat), 
  ##     -time matrix (time_mat2), 
  ##     -genetic matrix (MRCA_mat2), 
  ## compared to a reference:
  ##     -location matrix (geo_mat_ref), 
  ##     -time matrix (time_mat2_ref), 
  ##     -genetic matrix (MRCA_mat2_ref),
  ## for given bootstrap vectors x and y.
  
  geo_mat.tmp = geo_mat[x,x]
  time_mat.tmp2 = time_mat2[x,x]
  MRCA_mat.tmp2 = MRCA_mat2[x,x]
  geo_mat_ref.tmp = geo_mat_ref[y,y]
  time_mat2_ref.tmp2 = time_mat2_ref[y,y]
  MRCA_mat2_ref.tmp2 = MRCA_mat2_ref[y,y]
  
  tmp = MRCA_mat.tmp2 * time_mat.tmp2
  tmp[which(tmp == 0)] = NA
  
  tmp2 = time_mat.tmp2
  tmp2[which(tmp2 == 0)] = NA
  
  tmp_ref = MRCA_mat2_ref.tmp2 * time_mat2_ref.tmp2
  tmp_ref[which(tmp_ref == 0)] = NA
  
  tmp2_ref = time_mat2_ref.tmp2
  tmp2_ref[which(tmp2_ref == 0)] = NA
  
  a1 = cumsum(hist(tmp*geo_mat.tmp, breaks = c(0,Pmax,1E10), plot = F)$counts)
  a2 = cumsum(hist(tmp*geo_mat.tmp, breaks = c(0,Pmin,1E10), plot = F)$counts)
  a = a1 - a2

  b1 = cumsum(hist(tmp_ref * geo_mat_ref.tmp, breaks = c(0,Pmax,1E10), plot = F)$counts)
  b2 = cumsum(hist(tmp_ref * geo_mat_ref.tmp, breaks = c(0,Pmin,1E10), plot = F)$counts)
  b = b1 - b2

  c = sum(tmp2*geo_mat.tmp,na.rm=T)

  d = sum(tmp2_ref*geo_mat_ref.tmp,na.rm=T)
  rr.out = (a/c)/(b/d) 

  return(rr.out[-length(rr.out)])
}

#####################################################################
## Load relevant datasets and matrices
#####################################################################
## Metadata
data.seq = read.csv(args[1], sep = "\t")
## Time
time_mat = readRDS(args[2])
## Geography
geo_mat_country = readRDS(args[3])
geo_mat_continent = readRDS(args[4])
geo_mat_km_centroids = readRDS(args[5])
## Genetic distances
sim.mats <- readRDS(args[6])
nsim = dim(sim.mats)[3]

#####################################################################
## Parameters for the computation of the relative risks
#####################################################################
## Number of bootstrap event to perform, of each tree
nboot = 10

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
for (j in (1:n_steps)){
  print(j)
  # same country
  if (j == 1){
    for(ii in 1:nsim){
      print(ii)
      ## Choose isolates
      a = 1:nrow(data.seq)
      ## Choose location matrix: same country
      geo_mat = geo_mat_country[a,a]
      ## Time between isolates: max 2 years
      time_mat2 = time_mat[a,a]<=2
      ## Choose MRCA matrix
      MRCA_mat = sim.mats[,,ii]
      MRCA_mat2 = MRCA_mat[a,a]
      nseq = length(a)
      
      ## Reference
      ## Choose isolates: different countries >1000 km in the same continent
      ref = 1:nrow(data.seq)
      geo_mat_ref = (1-geo_mat_country[ref,ref])*(geo_mat_km_centroids[a,a]>1000)*(geo_mat_continent[ref,ref])
      geo_mat_ref[which(geo_mat_ref == 0)] = NA
      ## Time between isolates: max 2 years
      time_mat2_ref = time_mat[ref,ref]<=2
      ## Choose MRCA matrix
      MRCA_mat2_ref = MRCA_mat[ref,ref]
      nseq_ref = length(ref)
      
      ##Bootstrap to create the ci
      for (i in (1:nboot)){
        tmp = sample(nseq, nseq, replace = T)
        tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
        rr.out = ratio_bootstrap_dist_discrete_auto(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref)
        rr[((j*l-(l-1)):(j*l)),((ii-1)*nboot + i)] = rr.out
      }
    }
  }
  if (j == 2){ #different countries in the same continent, less than 1000km apart
    for(ii in 1:nsim){
      print(ii)
      ## Choose isolates
      a = 1:nrow(data.seq)
      ## Choose location matrix: different countries, <1000km apart, within the same continent
      geo_mat = (1-geo_mat_country[a,a])*(geo_mat_km_centroids[a,a]<=1000)*(geo_mat_continent[a,a])
      geo_mat[which(geo_mat == 0)] = NA
      ## Time between isolates: max 2 years
      time_mat2 = time_mat[a,a]<=2
      ## Choose MRCA matrix
      MRCA_mat2 = MRCA_mat[a,a]
      nseq = length(a)
      
      ## Reference
      ## Choose isolates: different countries >1000 km in the same continent
      ref = 1:nrow(data.seq)
      geo_mat_ref = (1-geo_mat_country[ref,ref])*(geo_mat_km_centroids[a,a]>1000)*(geo_mat_continent[ref,ref])
      geo_mat_ref[which(geo_mat_ref == 0)] = NA
      ## Time between isolates: max 2 years
      time_mat2_ref = time_mat[ref,ref]<=2
      MRCA_mat2_ref = MRCA_mat[ref,ref]
      nseq_ref = length(ref)
      
      ##Bootstrap to create the ci
      for (i in (1:nboot)){
        tmp = sample(nseq, nseq, replace = T)
        tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
        rr.out= ratio_bootstrap_dist_discrete_auto(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref)
        rr[((j*l-(l-1)):(j*l)),((ii-1)*nboot + i)] = rr.out
      }
    }
  }
  if (j == 3){ #different countries in the same continent, more than 1000km apart (reference)
    for(ii in 1:nsim){
      print(ii)
      ## Choose isolates
      a = 1:nrow(data.seq)
      ## Choose location matrix: different countries, >1000km apart, within the same continent
      geo_mat = (1-geo_mat_country[a,a])*(geo_mat_km_centroids[a,a]>1000)*(geo_mat_continent[a,a])
      geo_mat[which(geo_mat == 0)] = NA
      ## Time between isolates: max 2 years
      time_mat2 = time_mat[a,a]<=2
      ## Choose MRCA matrix
      MRCA_mat2 = MRCA_mat[a,a]
      
      ## Reference
      ## Choose isolates: different countries >1000 km in the same continent
      ref = 1:nrow(data.seq)
      geo_mat_ref = (1-geo_mat_country[ref,ref])*(geo_mat_km_centroids[a,a]>1000)*(geo_mat_continent[ref,ref])
      geo_mat_ref[which(geo_mat_ref == 0)] = NA
      ## Time between isolates: max 2 years
      time_mat2_ref = time_mat[ref,ref]<=2
      MRCA_mat2_ref = MRCA_mat[ref,ref]
      nseq_ref = length(ref)
      
      ##Bootstrap to create the ci
      for (i in (1:nboot)){
        tmp = sample(nseq, nseq, replace = T)
        tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
        rr.out= ratio_bootstrap_dist_discrete_auto(tmp, tmp, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref)
        rr[((j*l-(l-1)):(j*l)),((ii-1)*nboot + i)] = rr.out
      }
    }
  }
  if (j == 4){ #different continents
    for(ii in 1:nsim){
      print(ii)
      ## Choose isolates
      a = 1:nrow(data.seq)
      ## Choose location matrix: different continents
      geo_mat = 1-geo_mat_continent[a,a]
      ## Time between isolates: max 2 years
      time_mat2 = time_mat[a,a]<=2
      ## Choose MRCA matrix
      MRCA_mat2 = MRCA_mat[a,a]
      nseq = length(a)
      
      ## Reference
      ## Choose isolates: different countries >1000 km in the same continent
      ref = 1:nrow(data.seq)
      geo_mat_ref = (1-geo_mat_country[ref,ref])*(geo_mat_km_centroids[a,a]>1000)*(geo_mat_continent[ref,ref])
      geo_mat_ref[which(geo_mat_ref == 0)] = NA
      ## Time between isolates: max 2 years
      time_mat2_ref = time_mat[ref,ref]<=2
      ## Choose MRCA matrix
      MRCA_mat2_ref = MRCA_mat[ref,ref]
      nseq_ref = length(ref)
      
      ##Bootstrap to create the ci
      for (i in (1:nboot)){
        tmp = sample(nseq, nseq, replace = T)
        tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
        rr.out= ratio_bootstrap_dist_discrete_auto(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref)
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
## Write results
saveRDS(res, "risk_ratios.rds")
#####################################################################
