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
#' @param eval.index.grid see \code{DVDtest}
#' @param permat a result of \code{make.perm}
#' @param \dots partial arguments of \code{vdFun}
#' @param exclude an argument of \code{predict}
#' @param mc.cores a scalar, an argument of \code{mclapply}
#' @param dist.method see \code{DVDtest}
#' @return an array of distances
#' @author Meng Xu, Philip Reiss
#' @seealso \code{wass.perm}
#' @import parallel

wass_perm <-
function(vdFun, nperm, ydata1, ydata2, eval.index.grid, permat, ..., exclude,
                    mc.cores, dist.method){

  permarray <- array(dim = c(nperm, length(eval.index.grid), length(ydata1)))
  if (mc.cores == 1) {
    permtmp <- lapply(1:length(ydata1), wass.perm, vdFun = vdFun, dat1 = ydata1,
                      dat2 = ydata2, ..., permat = permat, .index = eval.index.grid, 
                      report.every = 10, exclude = exclude, dist.method = dist.method)
  for (k in 1:length(ydata1)) permarray[,,k] <- permtmp[[k]]
  return(permarray)
  }
  if (mc.cores != 1) {
    permtmp <- mclapply(1:length(ydata1), wass.perm, vdFun = vdFun, dat1 = ydata1, 
                        dat2 = ydata2, ..., permat=permat, .index=eval.index.grid, 
                        report.every = 10, exclude = exclude, mc.cores = mc.cores, 
                        dist.method = dist.method)
    for (k in 1:length(ydata1)) permarray[,,k] <- permtmp[[k]]
    return(permarray)
    }
  # permarray=array(dim=c(nperm,length(eval.index.grid),
  #                       length(ydata1)))
  # for (k in 1:length(ydata1)){
  #   permarray[,,k]<-wass.perm(k,vdFun=vdFun,dat1=ydata1[[k]], dat2=ydata2[[k]],
  #                             ...,
  #                             permat=permat.all, .index=eval.index.grid,
  #                             report.every=10, exclude=exclude)
  # }
  # return(permarray)
}
