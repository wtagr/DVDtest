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
    else if ("BCCG" %in% family) return(function(p) qBCCG(p, mu=mu, sigma=sigma, nu=nu))
    else if ("BCT" %in% family) return(function(p) qBCT(p, mu=mu, sigma=sigma, nu=nu, tau=tau))
    else if ("BCPE" %in% family) return(function(p) qBCPE(p, mu=mu, sigma=sigma, nu=nu, tau=tau))
  }
  if (any(dist.method == c('L2','L1','Hellinger'))) {
    if ("gaulss" %in% family) return(function(x) dnorm(x, mean=mu, sd=sigma))
    else if ("NO" %in% family) return(function(x) dnorm(x, mean=mu, sd=sigma))
    else if ("BCCG" %in% family) return(function(x) dBCCG(x, mu=mu, sigma=sigma, nu=nu))
    else if ("BCT" %in% family) return(function(x) dBCT(x, mu=mu, sigma=sigma, nu=nu, tau=tau))
    else if ("BCPE" %in% family) return(function(x) dBCPE(x, mu=mu, sigma=sigma, nu=nu, tau=tau))
  }
}
