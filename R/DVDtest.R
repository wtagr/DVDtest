#' Difference between Varying Distributions Test (DVDtest)
#' 
#' Testing the difference of two varying distributions.
#' 
#' This is the Details section
#' 
#' @param ydata1 a \code{data.frame} or a \code{list} of \code{data.frame}, containing 
#' at least 3 columns called '\code{.obs}', '\code{.index}' and '\code{.value}' which 
#' specify which curve the point belongs to (\code{.obs}) at which ('\code{.index}') 
#' it was observed and the observed value (\code{'.value'}). See details in the package 
#' \code{refund}. Other columns are available as well for modelling the varying 
#' distributions.
#' @param ydata2 same as \code{ydata1}. 
#' @param nperm a scalar, number of permutation
#' @param grid a vector, evaluation grids of \code{.index}
#' @param dist.method the distance measure to be used. This must be one of Wasserstein 
#' distance (\code{'wass'}), \code{'L2'} distance, \code{'L1'} distance and \code{'Hellinger'}. 
#' Defaults to \code{'wass'}.
#' @param mgcv.gam a logical variable, whether to apply \code{mgcv::gam} for eastimating 
#' distributions, whose parameters are a smooth function of a continuous variable. If 
#' \code{FALSE}, \code{gamlss.mx::gamlssNP} is adopted, which could cover a wider range of 
#' varing distributions.
#' @param \dots passed to arguments of \code{gam} or \code{gamlss}. If \code{mgcv.gam = TRUE}, 
#' \code{\dots} should include \code{formula}, \code{family} (=\code{gaulss()}) and 
#' other optional arguments in \code{mgcv::gam}. Otherwise, \code{...} passed to 
#' arguments inside of \code{gamlss::gamlss}.
#' @param exclude passed to \code{exclude} inside of \code{predict.gam} 
#' in case \code{mgcv.gam = TRUE}.
#' @param permadj a logical variable, whether to adjust the permutated data to cover 
#' the entire range, esp. in case of sparsity. Defaults to \code{FALSE}.
#' @param mc.cores passed to \code{mc.cores} inside of \code{mclapply} 
#' (not available on Windows unless \cr \code{mc.cores = 1}).The first element of \code{mc.cores} 
#' passed to dealing with the permutation process, the other passed to dealing with the long length 
#' of the data frame. Defaults to \code{c(1,1)}.
#' @param seeds set the seed for the permutation via \code{set.seed()}
#' @return 
#' \item{.index}{a vector, evaluation grids.}
#' \item{pval}{a vector or matrix of (corrected) p values.}
#' \item{vdparam}{a list of paramters of varying distributions.}
#' @note 
#' \itemize{
#' \item If \code{ydata1} and \code{ydata2} are \code{list}s of 
#' \code{data.frame}s, the lenghs of two lists must be the same.
#' 
#' \item If \code{mgcv.gam} is \code{TRUE}, \code{...} and \code{exclue} are \code{NULL} 
#' (default settings), then defaults to \cr
#' \code{formula <- list(.value ~ s(.index) + s(.obs, bs = "re"), ~ s(.index))} \cr
#' and \code{exclude <- "s(.obs)"}, repectively.
#' 
#' \item Normal distribution in mgcv::gam and BCCG, BCT and BCPE in 
#' gamlss::gamlss are currently supported by \code{DVDtest} for fitting a GAMLSS-type varying distributions. 
#' Please contact the maintainer if need further supporting.
#' }
#' @author Meng Xu, Philip Reiss
#' @references 
#' \itemize{
#' \item reiss-EMR18.pdf
#' \item Wood, S. N., Pya, N., & Säfken, B. (2016). Smoothing parameter and model selection for general 
#' smooth models. Journal of the American Statistical Association, 111(516), 1548-1563.
#' \item Rigby, R. A., & Stasinopoulos, D. M. (2005). Generalized additive models for location, scale 
#' and shape. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(3), 507-554.
#' }
#' 
#' @keywords permutation, pointwise
#' @export
#' @import mgcv gamlss gamlss.dist 
#' @examples
#' 
#' ## Data Generation ##
#'  p <- 6
#'  mu1 <- function(t) 0.2*(p-1)*sin(pi*t)+t+1
#'  mu2 <- function(t) -0.2*(p-1)*sin(pi*t)+t+1
#'  sig1 <- function(t) t+1
#'  sig2 <- sig1
#'  nperson <- 10
#'  fun1 <- function(t) rnorm(nperson, mu1(t), sig1(t))
#'  fun2 <- function(t) rnorm(nperson, mu2(t), sig2(t))
#'  tp <- seq(0,1,,10)
#'  data1 <- sapply(tp,fun1)
#'  data2 <- sapply(tp,fun2)
#'  
#'  library(reshape2)
#'  colnames(data2) <- colnames(data1) <- tp
#'  rownames(data2) <- 1:nperson+2*nperson
#'  dg1 <- melt(data1)
#'  dg2 <- melt(data2)
#'  colnames(dg1) <- colnames(dg2) <- c('.obs','.index','.value')
#'  dg1$.obs <- as.factor(dg1$.obs)
#'  dg2$.obs <- as.factor(dg2$.obs)
#'  # library(ggplot2)
#'  # ggplot() + geom_line(data = dg1, aes(x = .index,y = .value, col = factor(.obs)))
#'  # + geom_line(data = dg2, aes(x = .index, y = .value, col = factor(.obs)))
#' 
#'  ngrid <- 50
#'  ev.grid <- seq(0, 1, , ngrid)
#'  nperm. <- 50
#'  
#' ####Estimated with mgcv::gam
#'  simu.test1 <- DVDtest(dg1, dg2, nperm. ,ev.grid)
#' ####Estimated with gamlss::gamlss
#'  simu.test2 <- DVDtest(dg1, dg2, nperm.,ev.grid, formula = .value ~ cs(.index), 
#'                 sigma.formula = ~cs(.index), random = ~1|.obs, family = NO, mgcv.gam=F)
#' ####Plot
#'  simu.figs <- DVDplot(simu.test1)
#'  simu.figs$pfig
#'  simu.figs$kfig[[1]]



DVDtest <-
function(ydata1, ydata2, nperm, grid, dist.method = 'wass', mgcv.gam=TRUE,
               ..., exclude=NULL, permadj=FALSE, mc.cores = c(1,1), seeds = NULL){
  argmt <- list(...)
  if (!identical(sapply(ydata1,is.list), sapply(ydata2,is.list))) {
    stop('The types of ydata1 and ydata2 are not same')
  }
  if (!is.list(ydata1[[1]])) {
    ydata1<-list(ydata1)
    ydata2<-list(ydata2)
  }
  if (!all(c(c(".obs",".index",".value")%in%sapply(ydata1,names),
             c(".obs",".index",".value")%in%sapply(ydata2,names)))) {
    stop('The names of ydata1&2 need to include .index, .obs and .value at least')
  }
  if (length(ydata1) != length(ydata2)) stop('The lengths of ydata1 and ydata2 are not same')

  if (mgcv.gam & is.null(exclude)& !is.null(argmt[["formula"]])) stop('The arguments related mgcv::gam are missing, such as formula and exclude')

  if (mgcv.gam) vdFun<-mgcv::gam else vdFun<-gamlss::gamlss
  
  nroi <- length(ydata1)
  permat.all <- make.perms(ydata1[[1]], ydata2[[1]], nperm, grid, adj = permadj, seeds = seeds)
  reald <- get_realdist(vdFun = vdFun, ydata1 = ydata1, ydata2 = ydata2,
                        grid = grid,..., exclude = exclude, dist.method = dist.method,
                        mc.cores = mc.cores[2])
  realdists<-reald$realdists
  vdparam<-reald$vdparam
  permarray <- wass_perm(vdFun = vdFun, nperm = nperm, ydata1 = ydata1, ydata2 = ydata2,...,
                       grid = grid, permat = permat.all, exclude = exclude,
                       mc.cores = mc.cores[1], dist.method = dist.method)
  param.array <- get_params(nroi=nroi, nperm = nperm, permarray = permarray, grid = grid,
                            mc.cores = mc.cores[2])
  
  # Raw p-values (p.real, p.perm)
  p.mat <- get.pval(permarray = permarray, param.array = param.array, realdists = realdists,
                    nroi=nroi, grid = grid, nperm = nperm)
  return(list(.index = grid, pval = p.mat, vdparam = vdparam))
}
