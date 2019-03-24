#' Calculating the distances under the null hypothesis
#' 
#' 
#' @param vdFun a function, \code{gam} or \code{gamlss}, for fitting the varying distributions
#' @param ydata1 see \code{DVDtest}
#' @param ydata2 see \code{DVDtest}
#' @param ind.grid see \code{eval.index.grid} in \code{DVDtest}
#' @param \dots arguments of \code{vdFun}
#' @param excl an argument of \code{predict} for \code{gam}
#' @param mc.cores a scalar, an argument in \code{mclapply}
#' @param dist.method see \code{DVDtest}
#' @return a vector or matrix of the distances
#' @author Meng Xu, Philip Reiss
#' @seealso \code{DVDtest}


get_realdist <-
function(vdFun, ydata1, ydata2, ind.grid, ..., excl, mc.cores, dist.method){
  library(parallel)
  if (mc.cores == 1) {
    real.list <- sapply(1:length(ydata1), get.realdist, vdFun = vdFun,
         ydata1 = ydata1, ydata2 = ydata2,
         ind.grid = ind.grid, ...,excl = excl, dist.method = dist.method,
         simplify = TRUE)
    return(simplify2array(real.list))
  }
  if (mc.cores != 1) {
    real.list <- mclapply(1:length(ydata1), get.realdist, vdFun = vdFun,
                          ydata1 = ydata1, ydata2 = ydata2,
                          ind.grid = ind.grid, ..., excl = excl, mc.cores = mc.cores, 
                          dist.method = dist.method)
    realdists <- matrix(NA, length(ind.grid), length(ydata1))
    for (i in 1:length(ydata1)) realdists[,i] = real.list[[i]]
    return(realdists)
  }
}
