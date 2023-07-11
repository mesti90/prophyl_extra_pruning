#' Minimal Ancestor Deviation (MAD) rooting
#' 
#' @param unrooted_newick Unrooted tree string in newick format or a tree object
#' of class 'phylo'.
#' @param output_mode character; Amount of information to return. See Details
#' for more information.
#' @param cache logical; should computations be cached?
#' @param threads numeric; the number of CPU threads to use.
#' @param verbose logical; should verbose messages be printed to the console?
#' @return a list with the results containing one (\code{"newick"}), two
#' (\code{"stats"}) or six elements (\code{"full"})
#' @return if \code{cache} is \code{TRUE} also return a number of cache files:
#' \code{"distmat_t.rds"}, \code{"distmat_t2.rds"}, \code{"bad.rds"},
#' \code{"rho.rds"}.
#' @details If \code{output_mode} is \code{"newick"}(default), the function only
#' returns the rooted newick string. If \code{"stats"}, it also returns a
#' structure with the ambiguity index, clock cv, the minimum ancestor deviation
#' and the number of roots. If \code{"full"}, it also returns an unrooted tree
#' object, the index of the root branch, the branch ancestor deviations and a
#' rooted tree object.
#' @source Article: https://doi.org/10.1038/s41559-017-0193
#' @source Original code: https://www.mikrobio.uni-kiel.de/de/ag-dagan/ressourcen/mad2-2.zip
#' @examples
#' \dontrun {
#' # TODO simulate random tree
#' root_mad(unrooted_newick, output_mode)
#' }
#' @importFrom ape dist.nodes drop.tip is.binary is.rooted multi2di read.tree unroot 
root_mad <- function(unrooted_newick,
                     output_mode,
                     threads = 1,
                     cache = TRUE,
                     verbose = getOption("verbose")){
  output_mode <- match.arg(
    output_mode,
    choices = c("newick", "stats", "full")
  )
  t <- NA
  if(class(unrooted_newick)=="phylo"){ 
    t <- unrooted_newick
  }
  else{
    t <- ape::read.tree(text=unrooted_newick)
  }
  if(ape::is.rooted(t)){
    t <- ape::unroot(t)
  }
  #t$node.label<-NULL #To allow parsing when identical OTUs are present
  if(!ape::is.binary(t)){
    warning("Input tree is not binary! Internal multifurcations will be converted to branches of length zero and identical OTUs will be collapsed!")
    t <- ape::multi2di(t)
  }
  tf <- t$edge.length < 0
  if (any(tf)) {
    warning("Input tree contains negative branch lengths. They will be converted to zeros!")
    t$edge.length[tf]<-0
  }
  
  #### non-recursive alternative to collapse identical OTUs, if present
  collapsed_tree <- collapse_identical_tips(t)
  t <- collapsed_tree$tree
  dropped_tips <- collapsed_tree$dropped_tips
  #### End of non-recursive alternative
  
  notu <- length(t$tip.label)
  nbranch <- dim(t$edge)[1]
  npairs <- notu*(notu-1)/2
  nodeids <- 1:(nbranch+1)
  otuids <- 1:notu
  
  if (verbose) message("Calculating first distance matrix. ", appendLF = FALSE)
  if (cache && file.exists("distmat_t.rds")) {
    if (verbose) {
      message("Already calculated. Reading from cache. ", appendLF = FALSE)
    }
    dis <- readRDS("distmat_t.rds")
  } else {
    dis <- ape::dist.nodes(t) # phenetic distance. All nodes
    if (cache) {
      if (verbose) message("Saving to cache. ", appendLF = FALSE)
      saveRDS(dis, "distmat_t.rds")
    }
  }
  if (verbose) message("Done.")
  
  sdis <- dis[1:notu,1:notu] # phenetic distance. otus only

  #### Start recursion to collapse identical OTUs, if present.
  # ii<-which(sdis==0,arr.ind=TRUE)
  # k<-which(ii[,1]!=ii[,2])
  # if(length(k)){
  #   r<-ii[k[1],1]
  #   c<-ii[k[1],2]
  #   vv<-c(paste('@#',t$tip.label[r],'@#',sep=""),paste('(',t$tip.label[r],':0,',t$tip.label[c],':0)',sep=""))
  #   st<- ape::drop.tip(t,c) 
  #   st$tip.label[st$tip.label==t$tip.label[r]]<-vv[1]
  #   res<-root_mad(st,output_mode)
  #   if(is.list(res)){
  #     res[[1]]<-sub(vv[1],vv[2],res[[1]])
  #   }
  #   else{
  #     res<-sub(vv[1],vv[2],res)
  #   }
  #   return(res) #create the list 'res' to return the results 
  # }
  #### End of recursion


  
  t2 <- t
  t2$edge.length <- rep(1,nbranch)
  
  if (verbose) message("Calculating second distance matrix. ", appendLF = FALSE)
  if (cache && file.exists("distmat_t2.rds")) {
    if (verbose) {
      message("Already calculated. Reading from cache. ", appendLF = FALSE)
    }
    disbr <- readRDS("distmat_t2.rds")
  } else {
    disbr <- ape::dist.nodes(t2) # split distance. All nodes
    if (cache) {
      if (verbose) message("Saving to cache. ", appendLF = FALSE)
      saveRDS(dis, "distmat_t2.rds")
    }
  }
  if (verbose) message("Done.")
  
  sdisbr <- disbr[1:notu,1:notu] # split distance. otus only
  # NOTE
  # Before parallelisation, elements of this matrix were editted in a for loop.
  # Since parallelisation, elements of this matrix are only temporarily edited
  # within mclapply, but not globally. Is this a problem?
  i2p <- matrix(nrow = nbranch+1, ncol = notu)

  if (verbose) message("Calculating deviations. ", appendLF = FALSE)
  if (cache && file.exists("rho.rds") && file.exists("bad.rds")) {
    if (verbose) message("Already calculated. Reading from cache. ", appendLF = FALSE)
    rho <- readRDS("rho.rds")
    bad <- readRDS("bad.rds")
  } else {
    if (verbose) message("This may take a while. ", appendLF = FALSE)
    out <- parallel::mclapply(1:nbranch, function(br) {
      #collect the deviations associated with straddling otu pairs
      dij <- t$edge.length[br]
      if(dij==0){
        RHO<-NA
        BAD<-NA
        next
      }
      rbca <- numeric(npairs)
      i <- t$edge[br,1]
      j <- t$edge[br,2]
      sp <- dis[1:notu,i]<dis[1:notu,j] # otu split for 'br'
      dbc <- matrix(sdis[sp,!sp],nrow=sum(sp),ncol=sum(!sp))
      dbi <- replicate(dim(dbc)[2],dis[(1:notu)[sp],i]) 
      
      RHO <- sum((dbc-2*dbi)*dbc^-2)/(2*dij*sum(dbc^-2)) # optimized root node relative to 'i' node
      RHO <- min(max(0,RHO),1)
      dab <- dbi+(dij*RHO)
      ndab <- length(dab)
      rbca[1:ndab] <- as.vector(2*dab/dbc-1)
      # collect the remaining deviations (non-traversing otus)
      bcsp <- rbind(sp,!sp)
      ij <- c(i,j)
      counter <- ndab
      for (w in c(1,2)){
        if(sum(bcsp[w,])>=2){
          disbrw <- disbr[,ij[w]]
          pairids <- otuids[bcsp[w,]]
          for (z in pairids){
            i2p[,z] <- disbr[z,]+disbrw==disbrw[z]
          }
          for (z in 1:(length(pairids)-1)){
            p1 <- pairids[z]
            disp1 <- dis[p1,]
            pan <- nodeids[i2p[,p1]]
            for (y in (z+1):length(pairids)){
              p2 <- pairids[y]
              pan1 <- pan[i2p[pan,p2]]
              an <- pan1[which.max(disbrw[pan1])]
              counter <- counter+1
              rbca[counter] <- 2*disp1[an]/disp1[p2]-1
            }
          }
        }
      }
      if(length(rbca)!=npairs){
        stop("Unexpected number of pairs. Report this error to ftria@ifam.uni-kiel.de")
      }
      BAD <- sqrt(mean(rbca^2)) # branch ancestor deviation
      return(list(
        rho = RHO,
        bad = BAD
      ))
    }, mc.cores = threads)
    
    rho <- sapply(out, function(x) x$rho)
    bad <- sapply(out, function(x) x$bad)
    
    if (cache) {
      if (verbose) message("Saving to cache. ", appendLF = FALSE)
      saveRDS(rho, "rho.rds")
      saveRDS(bad, "bad.rds")
    }
  }
  if (verbose) message("Done.")
  
  # Select the branch with the minum ancestor deviation and calculate the root ambiguity index
  jj <- sort(bad,index.return = TRUE)
  tf<-bad==jj$x[1]
  tf[is.na(tf)]<-FALSE
  nroots <- sum(tf)
  if (nroots>1){
    warning("More than one possible root position. Multiple newick strings printed")
  }
  madr <- which(tf) # Index of the mad root branch(es)
  rai <- jj$x[1]/jj$x[2] # Root ambiguity index
  badr <- bad[tf] # Branch ancestor deviations value for the root(s)
  #Root the tree object, calculate the clock CV and retrieve the newick string
  rt <- vector(mode = "list",nroots) # Rooted tree object
  ccv <- vector(mode = "numeric",nroots) # Clock CV
  rooted_newick <- vector(mode = "character",nroots)
  
  for (i in 1:length(madr)){
    pp <- rho[madr[i]]*t$edge.length[madr[i]]
    nn <- t$edge[madr[i],]
    rt[[i]] <- phytools::reroot(t,nn[2],pos = pp)
    rooted_newick[i] <- ape::write.tree(rt[[i]])
    dd <- dis[1:notu,nn]
    sp <- dd[,1]<dd[,2]
    otu2root <- vector(mode="numeric",notu)
    otu2root[sp] <- dd[sp,1] + pp
    otu2root[!sp] <- dd[!sp,1] - pp
    ccv[i] <- 100*sd(otu2root)/mean(otu2root)
  }
  rooted_newick<-sub(')Root;',');',rooted_newick)
  # Output the result(s)
  if(missing(output_mode))
  {
    return(rooted_newick)
  }
  else{
    if(output_mode=='newick'){
      return(rooted_newick)
    }
    else if(output_mode=='stats'){ # Rooted newick and stats
      root_stats <- data.frame(ambiguity_index=rai,clock_cv=ccv,ancestor_deviation=badr,n_roots=nroots)
      return(list(rooted_newick,root_stats))
    }
    else if(output_mode=='full'){ #Rooted newick,stats, unrooted tree object, index of the branch root, ancestor deviations, rooted tree object
      root_stats <- data.frame(ambiguity_index=rai,clock_cv=ccv,ancestor_deviation=badr,n_roots=nroots)
      return(list(rooted_newick,root_stats,t,madr,bad,rt, dropped_tips))
    }
    else{
      return(rooted_newick)
    }
  }
}
