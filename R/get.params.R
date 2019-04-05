#' Fitting the generalized gamma (Pearson type III) distribution for each k
#' 
#' Fitting by the permuted-data distances
#' 
#' 
#' @param k a scalar, \code{k}th data.frame in \code{ydata1&2}
#' @param nperm a scalar, number of permutation
#' @param permarray an array, permuted-data distances from \code{wass_perm}
#' @param grid a vector, evaluation grids of \code{.index}
#' @return a list of mean, sigma, nu and the AIC
#' @author Philip Reiss, Meng Xu
#' @seealso \code{\link{get_params}}
#' @import gamlss gamlss.dist
#' @keywords internal
#' 
get.params <-
function(k, nperm, permarray, grid) {

  d1 <- na.omit(data.frame(dist = as.vector(t(permarray[,,k])), 
                           .index = rep(grid, nperm)))
  f1 <- quiet(gamlss(dist~pb(.index), sigma.formula = ~pb(.index), data = d1, family = GG))
  predd <- quiet(predictAll(f1, data = d1, newdata = data.frame(.index = grid)))
  list(mu = predd$mu[1:length(grid)], sigma = predd$sigma[1:length(grid)],
       nu = predd$nu[1:length(grid)], aic = f1$aic)
}
