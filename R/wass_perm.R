#' Calculating the distances via permuted data
#' 
#'
#' 
#'
#' 
#' @param vdFun a function, \code{gam} or \code{gamlss}, for fitting the varying distributions
#' @param nperm a scalar, number of permutation
#' @param ydata1 see \code{DVDtest}
#' @param ydata2 see \code{DVDtest}
#' @param grid see \code{\link{DVDtest}}
#' @param permat a result of \code{make.perm}
#' @param \dots partial arguments of \code{vdFun}
#' @param exclude an argument of \code{predict}
#' @param mc.cores a scalar, an argument of \code{mclapply}
#' @param dist.method see \code{DVDtest}
#' @return an array of distances
#' @author Meng Xu, Philip Reiss
#' @seealso \code{wass.perm}
#' @import parallel
#' @keywords internal

wass_perm <-
  function(vdFun, nperm, ydata1, ydata2, grid, permat, ..., exclude,
           mc.cores, dist.method){
    
    permarray <- array(dim = c(nperm, length(grid), length(ydata1)))
    
    permtmp <- mclapply(1:nperm, wass.perm, vdFun = vdFun, ydat1 = ydata1, 
                        ydat2 = ydata2, ..., permat=permat, grid=grid,
                        exclude = exclude, mc.cores = mc.cores, 
                        dist.method = dist.method)
    for (k in 1:nperm) permarray[k,,] <- permtmp[[k]]
    return(permarray)
  }