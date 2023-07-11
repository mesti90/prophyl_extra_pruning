collapse_identical_tips <- function(tree, seed = 0) {
  if (ape::is.rooted(tree)) {
    warning("Input tree is rooted. Unrooting.")
    tree <- ape::unroot(tree) 
  }
  if(!ape::is.binary(tree)){
    warning("Input tree is not binary. Converting to binary tree.")
    tree <- ape::multi2di(tree)
  }
  # Calculate cophenetic distances
  d <- ape::cophenetic.phylo(tree)
  
  set.seed(seed)
  
  drop_tbl <- data.frame()
  for (i in 1:nrow(d)) {
    index <- which(d[i,] == 0)
    if (length(index) > 1){
      index_keep <- sample(index, 1)
      index_drop <- index[-which(index == index_keep)]
      drop_tbl <- dplyr::bind_rows(
        drop_tbl,
        data.frame(
          drop = tree$tip.label[index_drop],
          keep = tree$tip.label[index_keep]
        )
      )
      d[index,index] <- NA
    }
  }
  diag(d) = 0
  
  tree <- ape::drop.tip(tree, tip = drop_tbl$drop)
  drop_tbl <- tibble::as_tibble(drop_tbl)
  
  return(list(
    tree = tree,
    dropped_tips = drop_tbl
  ))
}