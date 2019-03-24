#' Calculate the corrected P values
#' 
#' 
#' @param permarray an array, permuted-data distances from \code{wass_perm}
#' @param param.array an array, permuted-data distances from \code{get_params}
#' @param realdists a matrix or vector, the value from \code{get_realdist}
#' @param nroi a scalar, the length of \code{ydata1}
#' @param eval.index.grid a vector, evaluation grid of \code{.index}
#' @param nperm a scalar, number of permutation
#' @return a vector or matrix of p value
#' @author Philip Reiss, Meng Xu
#' @seealso \code{DVDtest}
#' 
get.pval <-
function(permarray, param.array, realdists, nroi, eval.index.grid, nperm){
  p.perm <- array(NA, dim = dim(permarray))
  p.real <- array(dim = dim(realdists))
  for (k in 1:length(eval.index.grid)) for (l in 1:nroi) {
    mu <- param.array[1,k,l]
    sigma <- param.array[2,k,l]
    nu <- param.array[3,k,l]
    p.real[k,l] <- pGG(realdists[k,l], mu, sigma, nu, lower.tail = FALSE)
    witch <- !is.na(permarray[,k,l])
    p.perm[witch,k,l] <- pGG(permarray[witch,k,l], mu, sigma,nu, lower.tail = FALSE)
  }
  
  # Adjusted p-values (pmat)
  minp <- apply(p.perm, 1, min, na.rm = TRUE)
  logicarray <- outer(minp, p.real, "<=")
  pmat <- (1 + apply(logicarray, 2:3, sum)) / (1 + nperm)
  return(pmat)
}
