#' Calculating the distances under the null hypothesis for each roi
#' 
#' 
#' @param k a scalar, \code{k}th data.frame in \code{ydata1&2}
#' @param vdFun a function, \code{gam} or \code{gamlss}, for fitting the varying distributions
#' @param ydata1 see \code{DVDtest}
#' @param ydata2 see \code{DVDtest}
#' @param ind.grid see \code{eval.index.grid} in \code{DVDtest}
#' @param \dots arguments of \code{vdFun}
#' @param excl an argument of \code{predict} for \code{gam}
#' @param dist.method see \code{DVDtest}
#' @return a vector of the distances
#' @author Philip Reiss, Meng Xu
#' @seealso \code{DVDtest}
#' @import mgcv
#' 
get.realdist <-
function(k, vdFun, ydata1, ydata2, ind.grid,
                       ..., excl, dist.method){

  argmt <- list(...)
  if (is.null(argmt[["formula"]])) {
    g1 <- vdFun(data = ydata1[[k]],
                formula = list(.value~s(.index)+s(.obs, bs="re"), ~s(.index)),
                family = gaulss())
    g2 <- vdFun(data=ydata2[[k]],
                formula = list(.value~s(.index)+s(.obs, bs="re"), ~s(.index)),
                family = gaulss())
    excl <- "s(.obs)"
  } else {
    g1 <- vdFun(data=ydata1[[k]],...)
    g2 <- vdFun(data=ydata2[[k]],...)
    }
  
  multiwass(g1, g2, newdata1=
              data.frame(.index=ind.grid, .obs=ydata1[[k]]$.obs[1]),
            newdata2 = data.frame(.index = ind.grid, .obs = ydata2[[k]]$.obs[1]),
            excl = excl, dist.method = dist.method)
}
