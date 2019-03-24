#' Distance functions
#' 
#' @param qfunc1 quantile or density functions from \code{params2qfunc}
#' @param qfunc2 quantile or density functions from \code{params2qfunc}
#' @param dist.method see \code{DVDtest}
#' @param \dots extra arguments in \code{integrate}
#' @return distance functions
#' @author Meng Xu, Philip Reiss
#' 
#' 
qfuncs2wass2 <-
function(qfunc1, qfunc2, dist.method = dist.method, ...) {
  if (dist.method=='wass'){
    diff.squared <- function(x) (qfunc1(x) - qfunc2(x)) ^2
    return(sqrt(integrate(diff.squared, 0, 1, ...)$value))
  }
  if (dist.method=='L2'){
    diff.squared <- function(x) (qfunc1(x) - qfunc2(x)) ^2
    return(sqrt(integrate(diff.squared, -Inf, Inf, ...)$value))
  }
  if (dist.method=='L1'){
    diff.squared <- function(x) abs(qfunc1(x) - qfunc2(x))
    return(sqrt(integrate(diff.squared, -Inf, Inf, ...)$value))
  }
  if (dist.method=='Hellinger'){
    diff.squared <- function(x) sqrt((qfunc1(x)*qfunc2(x)))
    return(1-sqrt(integrate(diff.squared, -Inf, Inf, ...)$value))
  }
}
