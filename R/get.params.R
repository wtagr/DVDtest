#' Fitting the generalized gamma (Pearson type III) distribution for each roi
#' 
#' Fitting by the permuted-data distances
#' 
#' 
#' @param roi a scalar, \code{roi}th data.frame in \code{ydata1&2}
#' @param nperm a scalar, number of permutation
#' @param permarray an array, permuted-data distances from \code{wass_perm}
#' @param eval.index.grid a vector, evaluation grid of \code{.index}
#' @return 
#' \item{mu}{a vector of mean}, \item{sigma}{a vector of sigma}, 
#' \item{nu}{a vector of nu}, \item{aic}{the Akaike information criterion in \code{gamlss}}
#' @author Philip Reiss, Meng Xu
#' @seealso \code{\link{get_params}}
#' 
#' 
get.params <-
function(roi, nperm, permarray, eval.index.grid) {
  library(gamlss)
  d1 <- na.omit(data.frame(dist = as.vector(t(permarray[,,roi])), 
                           .index = rep(eval.index.grid, nperm)))
  f1 <- gamlss(dist~pb(.index), sigma.fo=~pb(.index), data=d1, family=GG)
  predd <- predictAll(f1, data = d1, newdata = data.frame(.index = eval.index.grid))
  list(mu = predd$mu[1:length(eval.index.grid)], sigma = predd$sigma[1:length(eval.index.grid)],
       nu = predd$nu[1:length(eval.index.grid)], aic = f1$aic)
}
