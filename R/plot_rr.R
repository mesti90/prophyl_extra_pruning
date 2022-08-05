#' Plot Relative Risks
#' 
#' This function plots the results of a risk analysis carried out within the
#' analysis pipeline.
#' @param dir character; the directory from which the pipeline was launched
#' @return The function returns a plot.
#' @export
plot_rr <- function(dir = "."){
  
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
  
  res <- readRDS(paste0(dir, "/results/calculate_risk_ratios/risk_ratios.rds"))
  int = res$int
  l <- length(int) -1
  nboot = res$nboot
  nsim = res$nsim
  rr <- res$rr
  
  ## Plot
  rr[which(is.infinite(rr))] = NA
  par(mfrow = c(l,1), mar = c(2, 4, 1, 0.5), oma = c(7, 2, 0, 0))
  n_steps = 4
  for (i in (1:l)){
    boot.ci = apply(rr[c(0:(n_steps-1)*(l)+i),(1:(nboot*nsim))], 1, quantile, probs = c(0.025,0.975), na.rm = T)
    boot.ci.m2 = apply(rr[c(0:(n_steps-1)*(l)+i),(1:(nboot*nsim))], 1, quantile, probs = c(0.5), na.rm = T)
    print(boot.ci)
    print(boot.ci.m2)
    
    boot.ci[which(boot.ci > 1E7)] = 1E7
    boot.ci[which(boot.ci < 0.01)] = 0.01
    boot.ci.m2[which(boot.ci.m2 < 0.01)] = 0.01
    
    plot(1:(n_steps), boot.ci.m2, type="p", pch=16,cex=1.4, xlab="", xlim = c(0.5,n_steps+0.5),
         ylim = c(0.01,500), log='y',
         ylab="", main = paste0(int[i], " <MRCA< ", int[i+1], " years"), xaxt="n", cex.axis= 1.1, yaxt="n")
    axis(2, at = c(0.01, 0.1, 1, 10, 100, 1000), label =  c('<0.01', 0.1, 1, 10, 100 ,1000), las = 2, cex.axis= 1.1)
    if(i == l){
      axis(1, at = 1:(n_steps), labels = c("Wihtin countries",
                                           "Between countries \n <1000km",
                                           "Between countries \n >1000km (ref)",
                                           "Between continents"), cex.axis=1, las =2)
    }
    else{
      axis(1, at = 1:(n_steps), labels = c("", "", "", ""), cex.axis=1, las =2)
    }
    abline(h = 1, col = 'red')
    arrows(1:(n_steps), boot.ci[1,], 1:(n_steps), boot.ci[2,], length=0.05, angle=90, code=3)
  }
  mtext(paste0("Relative risk"), outer = T, cex = 1.3, side = 2)
}
