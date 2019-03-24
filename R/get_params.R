#' Fitting the generalized gamma (Pearson type III) distribution
#' 
#'
#' 
#'
#' 
#' @param nroi a scalar, the length of \code{ydata1} or \code{ydata2}
#' @param nperm a scalar, number of permutation
#' @param permarray an array, permuted-data distances from \code{wass_perm}
#' @param eval.index.grid a vector, evaluation grid of \code{.index}
#' @return an array, \code{param.array}
#' @note NULL
#' @author Meng Xu, Philip Reiss
#' @seealso \code{DVDtest}
#' 
#' 
get_params <-
function(nroi, nperm, permarray, eval.index.grid){
  param.array<- array(dim = c(3,length(eval.index.grid), nroi))
  require(gamlss)
  predlist <- lapply(1:nroi, get.params, nperm = nperm, permarray = permarray,
                   eval.index.grid = eval.index.grid)
  for (roi in 1:nroi) {
    param.array[1,,roi] <- predlist[[roi]]$mu
    param.array[2,,roi] <- predlist[[roi]]$sigma
    param.array[3,,roi] <- predlist[[roi]]$nu
  }
  return(param.array)
}
