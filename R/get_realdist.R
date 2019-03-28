#' Calculating the distances under the null hypothesis
#' 
#' 
#' @param vdFun a function, \code{gam} or \code{gamlss}, for fitting the varying distributions
#' @param ydata1 see \code{DVDtest}
#' @param ydata2 see \code{DVDtest}
#' @param grid see in \code{DVDtest}
#' @param \dots arguments of \code{vdFun}
#' @param excl an argument of \code{predict} for \code{gam}
#' @param mc.cores a scalar, an argument in \code{mclapply}
#' @param dist.method see \code{DVDtest}
#' @return a vector or matrix of the distances and list of parameters of varying 
#' distributions.
#' @author Meng Xu, Philip Reiss
#' @seealso \code{\link{DVDtest}}
#' @import parallel
#' @keywords internal


get_realdist <-
  function(vdFun, ydata1, ydata2, grid, ..., excl, mc.cores, dist.method){
    vd_param <- list()
    realdists <- matrix(NA, length(grid), length(ydata1))
    real.list <- mclapply(1:length(ydata1), get.realdist, vdFun = vdFun,
                          ydata1 = ydata1, ydata2 = ydata2,
                          grid = grid, ..., excl = excl, mc.cores = mc.cores, 
                          dist.method = dist.method)
    for (i in 1:length(ydata1)) {
      realdists[,i] <- real.list[[i]]$rlist
      vd_param[[i]] <- real.list[[i]]$vd.param
    }
    return(list(realdists=realdists,vdparam=vd_param))
  }
