#' Calculating the distances under the null hypothesis for each k
#' 
#' 
#' @param k a scalar, \code{k}th data.frame in \code{ydata1&2}
#' @param ydata1 see \code{DVDtest}
#' @param ydata2 see \code{DVDtest}
#' @param grid see \code{DVDtest}
#' @param \dots arguments of \code{vdFun}
#' @param excl an argument of \code{predict} for \code{gam}
#' @param dist.method see \code{DVDtest}
#' @return a vector of the distances and parameters of varying distributions
#' @author Philip Reiss, Meng Xu
#' @seealso \code{\link{DVDtest}}
#' @import mgcv gamlss
#' @keywords internal
#' 
get.realdist <- function(k, ydata1, ydata2, grid, ..., exclude, dist.method, mgcv.gam) {
  
  argmt <- list(...)
  ydat1 <- ydata1[[k]]
  ydat2 <- ydata2[[k]]

  if (is.null(argmt[["formula"]])) {
    g1 <- quiet(gam(data = ydat1, formula = list(.value ~ s(.index) + 
                                                         s(.obs, bs = "re"), ~s(.index)), family = gaulss()))
    g2 <- quiet(gam(data = ydat2, formula = list(.value ~ s(.index) + 
                                                         s(.obs, bs = "re"), ~s(.index)), family = gaulss()))
    exclude <- "s(.obs)"
  } else if (mgcv.gam) {
    g1 <- quiet(gam(data = ydat1, formula = argmt[["formula"]], 
                    family = gaulss()))
    g2 <- quiet(gam(data = ydat2, formula = argmt[["formula"]], 
                    family = gaulss()))
    exclude <- "s(.obs)"
  } else {
    g1 <- quiet(gamlss(data = ydat1, ...))
    g2 <- quiet(gamlss(data = ydat2, ...))
  }
  
  rlist <- multiwass(g1, g2, newdata1 = data.frame(.index = grid, .obs = ydat1$.obs[1]), 
                     newdata2 = data.frame(.index = grid, .obs = ydat2$.obs[1]), 
                     exclude = exclude, dist.method = dist.method, dt1 = ydat1, 
                     dt2 = ydat2)
  return(list(rlist = rlist$wvec, vd.param = rlist$pred))
}
