#' Getting the quantile/density function via parameters
#' 
#' @param params a vector, parameters for certain distributions
#' @param family a specific distribution from \code{gam} or \code{gamlss}
#' @param dist.method see \code{DVDtest}
#' @return quantile or density functions
#' @author Meng Xu, Philip Reiss
#' @import gamlss gamlss.dist
#' @keywords internal
params2qfunc <-
function(params, family, dist.method) {
  mu <- params$mu
  sigma <- params$sigma
  nu <- params$nu
  tau <- params$tau
  if (dist.method == 'wass'){
    if ("gaulss" %in% family) return(function(p) qnorm(p, mean=mu, sd=sigma))
    else if ("NO" %in% family) return(function(p) qnorm(p, mean=mu, sd=sigma))
    else {
      qfun. <- match.fun(family[1])
      qfun <- paste0("q",qfun.)
      return(function(p) qfun(p, mu=mu, sigma=sigma, nu=nu, tau=tau))
    }
  }

  if (any(dist.method == c('L2','L1','Hellinger'))) {
    if ("gaulss" %in% family) return(function(x) dnorm(x, mean=mu, sd=sigma))
    else if ("NO" %in% family) return(function(x) dnorm(x, mean=mu, sd=sigma))
    else {
      pfun. <- match.fun(family[1])
      pfun <- paste0("d",pfun.)
      return(function(x) pfun(x, mu=mu, sigma=sigma, nu=nu, tau=tau))
    }
  }
  }