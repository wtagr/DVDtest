#' Calculating the distances via permuted data for each \code{k}
#' 
#' 
#' @param k a scalar, \code{k}th \code{data.frame} of \code{ydata1&2}
#' @param vdFun a function, \code{gam} or \code{gamlss}, for fitting the varying distributions
#' @param dat1 \code{k}th \code{data.frame} of \code{ydata1}
#' @param dat2 \code{k}th \code{data.frame} of \code{ydata2}
#' @param \dots arguments of \code{vdFun}
#' @param permat a result of \code{make.perm}
#' @param .index see \code{eval.index.grid} in \code{DVDtest}
#' @param report.every a scalar, reporting the number permutation
#' @param exclude an argument of \code{predict}  
#' @param dist.method see \code{DVDtest}
#' @return a matrix of permuted-data distances
#' @author Philip Reiss, Meng Xu
#' @seealso \code{wass_perm}
#' @import mgcv
#' 

wass.perm <-
function(k, vdFun, dat1, dat2, ..., permat, .index, 
                      report.every = 10, exclude, dist.method) {

  dat1 <- dat1[[k]]
  dat2 <- dat2[[k]]
  n1 <- length(unique(dat1$.obs))
  n2 <- length(unique(dat2$.obs))
  nerror <- 0
  nperm <- nrow(permat)
  permdists <- matrix(NA, nperm, length(.index))
  argmt <- list(...)
  for (i in 1:nperm) {
    if (!i%%report.every) cat("*** Perm", i, "***\n")
    bothdat <- rbind(dat1, dat2)
    which1 <- permat[i, ]
    perm1 <- bothdat[bothdat$.obs %in% unique(bothdat$.obs)[which1],]
    perm2 <- bothdat[!bothdat$.obs %in% unique(bothdat$.obs)[which1],]
    
    if (is.null(argmt[["formula"]])) {
      g1.p <- vdFun(data = perm1,
                  formula = list(.value~s(.index)+s(.obs, bs="re"), ~s(.index)),
                  family = gaulss())
      g2.p <- vdFun(data = perm2,
                  formula = list(.value~s(.index)+s(.obs, bs="re"), ~s(.index)),
                  family = gaulss())
      excl <- "s(.obs)"
    } else {
      g1.p <- try(vdFun(data = perm1, ...))
      g2.p <- try(vdFun(data = perm2, ...))
    }
    
    if (inherits(g1.p, "try-error") || inherits(g2.p, "try-error")) nerror <- nerror + 1
    else permdists[i, ] <- multiwass(g1.p, g2.p, newdata1=data.frame(.index=.index, .obs=perm1$.obs[1]), newdata2=data.frame(.index=.index, .obs=perm2$.obs[1]),
                                     dist.method = dist.method, exclude=exclude) 
  }
  # list(.index=.index, nperm=nperm, nerror=nerror, permdists=permdists)
  return(permdists)
}
