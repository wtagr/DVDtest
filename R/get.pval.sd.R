#' Calculate the corrected P values with step-down process
#' 
#' 
#' @param permarray an array, permuted-data distances from \code{wass_perm}
#' @param param.array an array, permuted-data distances from \code{get_params}
#' @param realdists a matrix or vector, the value from \code{get_realdist}
#' @param nroi a scalar, the length of \code{ydata1}
#' @param grid a vector, evaluation grids of \code{.index}
#' @param nperm a scalar, number of permutation
#' @return a vector or matrix of p value
#' @author Philip Reiss, Meng Xu
#' @seealso \code{\link{DVDtest}}
#' @import gamlss.dist
#' @keywords internal
#' 
get.pval.sd <- function(permarray, param.array, realdists, nroi, grid, 
                        nperm) {
  p.perm <- array(NA, dim = dim(permarray))
  p.real <- array(dim = dim(realdists))
  for (k in 1:length(grid)) for (l in 1:nroi) {
    mu <- param.array[1, k, l]
    sigma <- param.array[2, k, l]
    # nu <- param.array[3, k, l]
    p.real[k, l] <- pGA(realdists[k, l], mu, sigma, lower.tail = FALSE)
    witch <- !is.na(permarray[, k, l])
    p.perm[witch, k, l] <- pGA(permarray[witch, k, l], mu, sigma, lower.tail = FALSE)
  }
  
  p.r <- as.vector(matrix(p.real))
  p.p <- matrix(NA, nperm, length(p.r))
  for (i in 1:nperm) p.p[i, ] <- matrix(p.perm[i, , ])
  
  p.rrk <- order(p.r)
  p.rstar <- sort(p.r)
  # q.pstar.t<-t(apply(p.p,1,sort))
  q.pstar <- p.p
  for (b in 1:nperm) for (i in 1:length(p.r)) q.pstar[b, i] = min(p.p[b, 
                                                                      i:length(p.r)], na.rm = T)
  
  qmat <- q.pstar
  for (b in 1:nperm) for (i in (length(p.r) - 1):1) qmat[b, i] = min(qmat[b, 
                                                                          i + 1], q.pstar[b, i], na.rm = T)
  
  p.tilda.star = p.r
  for (i in 1:length(p.r)) p.tilda.star[i] <- sum(qmat[, i] <= p.rstar[i], 
                                                  na.rm = T)/(nperm + 1)
  # monotonicity
  for (i in 2:length(p.r)) p.tilda.star[i] = max(p.tilda.star[i - 1], 
                                                 p.tilda.star[i], na.rm = T)
  
  p.tilda = p.tilda.star[order(p.rrk)]
  pmat.sd <- matrix(p.tilda, length(grid))
  return(pmat.sd)
}
