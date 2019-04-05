#' Calculating the distances between two gam/gamlss objects
#' 
#' 
#' 
#' @param obj1 a gam/gamlss object
#' @param obj2 another gam/gamlss object
#' @param newdata1 related evaluation grids
#' @param newdata2 related evaluation grids
#' @param dt1 raw data
#' @param dt2 raw data
#' @param dist.method see \code{DVDtest}
#' @param \dots partial arguments in \code{predict}
#' @return a vector, distances
#' @author Philip Reiss, Meng Xu
#' @keywords internal
multiwass <-
function(obj1, obj2, newdata1, newdata2, dist.method, dt1, dt2, ...) {
  npts <- NROW(newdata1)
  wvec <- c()
  if ("gamlss" %in% c(class(obj1))) pred1 <- quiet(predictAll(obj1, newdata1, data = dt1))
  else if ("gaulss" %in% obj1$family) {
    tmp <- quiet(predict(obj1, newdata1, type = "response", ...))
    pred1 <- list(mu = tmp[,1], sigma = 1/tmp[,2])
  }
  if ("gamlss" %in% c(class(obj2))) pred2 <- quiet(predictAll(obj2, newdata2, data = dt2))
  else if ("gaulss" %in% obj2$family) {
    tmp <- quiet(predict(obj2, newdata2, type = "response", ...))
    pred2 <- list(mu = tmp[,1], sigma = 1/tmp[,2])
  }
  for (k in 1:npts) {
    params1 <- list(mu = pred1$mu[k], sigma = pred1$sigma[k], nu = pred1$nu[k], tau = pred1$tau[k])
    params2 <- list(mu = pred2$mu[k], sigma = pred2$sigma[k], nu = pred2$nu[k], tau = pred2$tau[k])
    qf1 <- params2qfunc(params1, obj1$family, dist.method = dist.method)
    qf2 <- params2qfunc(params2, obj2$family, dist.method = dist.method)
    wvec[k] <- qfuncs2wass2(qf1, qf2, dist.method = dist.method)
  }
  if (NCOL(newdata1) == 1) names(wvec) = newdata1[,1]
  return(list(wvec = wvec, pred = list(pred1 = pred1, pred2 = pred2)))  # can we make this a *named* vector?
}
