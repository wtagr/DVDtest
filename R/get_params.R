#' Fitting the Gamma distribution (Pearson type III) 
#' 
#'
#' 
#'
#' 
#' @param nroi a scalar, the length of \code{ydata1} or \code{ydata2}
#' @param nperm a scalar, number of permutation
#' @param permarray an array, permuted-data distances from \code{wass_perm}
#' @param grid a vector, evaluation grid of \code{.index}
#' @param mc.cores a number of the cores using for multi-process
#' @return an array, \code{param.array}
#' @note NULL
#' @author Meng Xu, Philip Reiss
#' @seealso \code{DVDtest}
#' @import gamlss
#' @keywords internal
get_params <-
function(nroi, nperm, permarray, grid, mc.cores){
  param.array<- array(dim = c(2, length(grid), nroi))

  predlist <- mclapply(1:nroi, get.params, nperm = nperm, permarray = permarray,
                   grid = grid, mc.cores = mc.cores)
  for (roi in 1:nroi) {
    param.array[1,,roi] <- predlist[[roi]]$mu
    param.array[2,,roi] <- predlist[[roi]]$sigma
    # param.array[3,,roi] <- predlist[[roi]]$nu
  }
  return(param.array)
}
