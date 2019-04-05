#' Calculating the distances via permuted data for each k
#' 
#' 
#' @param k a scalar, \code{k}th \code{data.frame} of \code{ydata1&2}
#' @param vdFun a function, \code{gam} or \code{gamlss}, for fitting the varying distributions
#' @param dat1 \code{k}th \code{data.frame} of \code{ydata1}
#' @param dat2 \code{k}th \code{data.frame} of \code{ydata2}
#' @param \dots arguments of \code{vdFun}
#' @param permat a result of \code{make.perm}
#' @param grid see \code{grid} in \code{DVDtest}
#' @param exclude an argument of \code{predict}  
#' @param dist.method see \code{DVDtest}
#' @return a matrix of permuted-data distances
#' @author Philip Reiss, Meng Xu
#' @seealso \code{wass_perm}
#' @import mgcv
#' @keywords internal

wass.perm <-
  function(np, vdFun, ydat1, ydat2, ..., permat, grid, exclude, dist.method) {
    which1 <- permat[np, ]
    nroi <- length(ydat1)
    permdists<-matrix(NA,length(grid),nroi)
    
    for (k in 1:nroi){
      dat1 <- ydat1[[k]]
      dat2 <- ydat2[[k]]
      n1 <- length(unique(dat1$.obs))
      n2 <- length(unique(dat2$.obs))
      nerror <- 0
      nperm <- nrow(permat)
      
      argmt <- list(...)
      bothdat <- rbind(dat1, dat2)
      
      perm1 <- bothdat[bothdat$.obs %in% unique(bothdat$.obs)[which1],]
      perm2 <- bothdat[!bothdat$.obs %in% unique(bothdat$.obs)[which1],]
      if (is.null(argmt[["formula"]])) {
        g1.p <- try(vdFun(data = perm1,
                      formula = list(.value~s(.index)+s(.obs, bs="re"), ~s(.index)),
                      family = gaulss()))
        g2.p <- try(vdFun(data = perm2,
                      formula = list(.value~s(.index)+s(.obs, bs="re"), ~s(.index)),
                      family = gaulss()))
        exclude <- "s(.obs)"
      } else {
        g1.p <- try(vdFun(data = perm1, ...))
        g2.p <- try(vdFun(data = perm2, ...))
      }
      
      if (inherits(g1.p, "try-error") || inherits(g2.p, "try-error")) nerror <- nerror + 1
      else permdists[,k] <- multiwass(g1.p, g2.p, newdata1=data.frame(.index=grid, .obs=perm1$.obs[1]), newdata2=data.frame(.index=grid, .obs=perm2$.obs[1]),
                                      dist.method = dist.method, exclude=exclude)$wvec 
    }
    # list(.index=.index, nperm=nperm, nerror=nerror, permdists=permdists)
    return(permdists)
  }
