#' Calculating the distances under the null hypothesis for each k
#' 
#' 
#' @param k a scalar, \code{k}th data.frame in \code{ydata1&2}
#' @param vdFun a function, \code{gam} or \code{gamlss}, for fitting the varying distributions
#' @param ydata1 see \code{DVDtest}
#' @param ydata2 see \code{DVDtest}
#' @param grid see \code{DVDtest}
#' @param \dots arguments of \code{vdFun}
#' @param excl an argument of \code{predict} for \code{gam}
#' @param dist.method see \code{DVDtest}
#' @return a vector of the distances and parameters of varying distributions
#' @author Philip Reiss, Meng Xu
#' @seealso \code{\link{DVDtest}}
#' @import mgcv
#' @keywords internal
#' 
get.realdist <-
function(k, vdFun, ydata1, ydata2, grid, ..., exclude, dist.method){

  argmt <- list(...)
  if (is.null(argmt[["formula"]])) {
    g1 <- quiet(vdFun(data = ydata1[[k]],
                formula = list(.value~s(.index)+s(.obs, bs="re"), ~s(.index)),
                family = gaulss()))
    g2 <- quiet(vdFun(data=ydata2[[k]],
                formula = list(.value~s(.index)+s(.obs, bs="re"), ~s(.index)),
                family = gaulss()))
    exclude <- "s(.obs)"
  } else {
    g1 <- quiet(vdFun(data=ydata1[[k]],...))
    g2 <- quiet(vdFun(data=ydata2[[k]],...))
    }
  
  rlist <- multiwass(g1, g2, newdata1 = data.frame(.index = grid, .obs = ydata1[[k]]$.obs[1]),
            newdata2 = data.frame(.index = grid, .obs = ydata2[[k]]$.obs[1]),
            exclude = exclude, dist.method = dist.method, dt1 = ydata1[[k]], dt2 = ydata2[[k]])
  return(list(rlist = rlist$wvec, vd.param = rlist$pred))
}
