#' Calculating the distances under the null hypothesis
#' 
#' 
#' @param ydata1 see \code{DVDtest}
#' @param ydata2 see \code{DVDtest}
#' @param grid see in \code{DVDtest}
#' @param \dots arguments of \code{vdFun}
#' @param excl an argument of \code{predict} for \code{gam}
#' @param dist.method see \code{DVDtest}
#' @param mc.cores see the 2nd element of \code{mc.cores} in \code{DVDtest}
#' @return a vector or matrix of the distances and list of parameters of varying 
#' distributions.
#' @author Meng Xu, Philip Reiss
#' @seealso \code{\link{DVDtest}}
#' @import parallel
#' @keywords internal


get_realdist <- function(ydata1, ydata2, grid, ..., exclude, mc.cores, 
                         dist.method, mgcv.gam) {
  vd_param <- list()
  realdists <- matrix(NA, length(grid), length(ydata1))
  real.list <- mclapply(1:length(ydata1), get.realdist, ydata1 = ydata1, 
                        ydata2 = ydata2, grid = grid, ..., exclude = exclude, dist.method = dist.method, 
                        mc.cores = mc.cores, mgcv.gam = mgcv.gam)
  for (i in 1:length(ydata1)) {
    realdists[, i] <- real.list[[i]]$rlist
    vd_param[[i]] <- real.list[[i]]$vd.param
  }
  return(list(realdists = realdists, vdparam = vd_param))
}

