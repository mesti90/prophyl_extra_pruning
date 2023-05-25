#' Calculate relative risks for a geographic matrix
#' 
#' @param x integer; vector of tree tip indices. If tree is not bootstrapped, 
#' a sequence from 1 until the number of tips within the tree. If tree is
#' bootstrapped, a bootstrap sample of this sequence (equal length with
#' replacement).
#' @param y integer; vector of tree tip indices. If tree is not bootstrapped, 
#' a sequence from 1 until the number of tips within the tree. If tree is
#' bootstrapped, a bootstrap sample of this sequence (equal length with
#' replacement).
#' @param geo_mat matrix; a symmetric geographic matrix, rows and columns are
#' about tips of the tree, cells contain geographical information, e.g. whether
#' the isolates for the two tips are from the same country.
#' @param time_mat2 matrix; a symmetric logical matrix, rows and columns are
#' about tips of the tree, cells tell whether the difference between the
#' sampling dates of two isolates are below a certain threshold, e.g. 2 years.
#' Pairs which pass receive value 1. pairs that don't receive value 0.
#' @param MRCA_mat2 matrix; a symmetric temporal matrix, rows and columns are
#' about tips of the tree, cells contain the cophenetic distance of tips, in
#' years.
#' @param geo_mat_ref matrix; symmetric logical matrix, rows and columns are
#' about tips of a tree, cells indicate whether pair is a geographic reference
#' pair.
#' @param time_mat2_ref matrix; symmetric logical matrix, rows and columns are
#' about tips of a tree, cells indicate whether the difference between the
#' sampling dates of two isolates are below a certain threshold, e.g. 2 years.
#' Pairs which pass receive value 1. pairs that don't receive value 0.
#' @param MRCA_mat2_ref matrix; a symmetric temporal matrix, rows and columns are
#' about tips of the tree, cells contain the cophenetic distance of tips, in
#' years.
#' @param Pmax numeric; vector of max temporal windows
#' @param Pmin numeric; vector of min temporal windows
#' @note within the pipeline application it seems \code{time_mat2} and
#' \code{time_mat2_ref} are the same, and similarly, \code{MRCA_mat2} and
#' \code{MRCA_mat2_ref} are the same. TODO: check and simplify if possible.
ratio_bootstrap_dist_discrete_auto <- function(
    x,
    y,
    geo_mat,
    time_mat2,
    MRCA_mat2,
    geo_mat_ref,
    time_mat2_ref,
    MRCA_mat2_ref,
    Pmax,
    Pmin) {
  
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