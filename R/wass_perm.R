#' Calculating the distances via permuted data
#' 
#'
#' 
#'
#' 
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


wass_perm <- function(nperm, ydata1, ydata2, grid, permat, ..., exclude, 
                      mc.cores, dist.method, mgcv.gam) {
  
  permarray <- array(dim = c(nperm, length(grid), length(ydata1)))
  
  permtmp <- mclapply(1:nperm, wass.perm, ydat1 = ydata1, 
                      ydat2 = ydata2, ..., permat = permat, grid = grid, exclude = exclude, 
                      dist.method = dist.method, mgcv.gam = mgcv.gam, 
                      mc.cores = mc.cores)
  for (k in 1:nperm) permarray[k, , ] <- permtmp[[k]]
  return(permarray)
}