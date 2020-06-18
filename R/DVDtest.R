#' Difference between Varying Distributions Test (DVDtest)
#' 
#' Testing the difference of two varying distributions.
#' 
#' This is the Details section
#' 
#' @param ydata1 a \code{data.frame} or a \code{list} of \code{data.frame}, containing 
#' at least 3 columns called \code{.obs}, \code{.index} and \code{.value} which 
#' specify which curve the point belongs to (\code{.obs}) at which (\code{.index}) 
#' it was observed and the observed value (\code{.value}). See details in the package 
#' \code{refund}. Other columns are available as well for modelling the varying 
#' distributions.
#' @param ydata2 same as \code{ydata1}. 
#' @param nperm a scalar, number of permutation
#' @param grid a vector, evaluation grids of \code{.index}
#' @param dist.method the distance measure to be used. This must be one of Wasserstein 
#' distance (\code{'wass'}), \code{'L2'} distance, \code{'L1'} distance and \code{'Hellinger'}. 
#' Defaults to \code{'wass'}.
#' @param mgcv.gam a logical variable, whether to apply \code{\link[mgcv]{gam}} for eastimating 
#' distributions, whose parameters are a smooth function of a continuous variable. If 
#' \code{FALSE}, \code{\link[gamlss]{gamlss}} is adopted, which could cover a wider range of 
#' varying distributions.
#' @param \dots passed to arguments of \code{\link[mgcv]{gam}} or \code{\link[gamlss]{gamlss}}. If \code{mgcv.gam = TRUE}, 
#' \code{\dots} should include \code{formula}, \code{family} and 
#' other optional arguments in \code{mgcv::gam}. Otherwise, \code{...} passed to 
#' arguments inside of \code{\link[gamlss]{gamlss}}. See Examples for details.
#' @param exclude passed to \code{exclude} inside of \code{predict.gam} 
#' in case \code{mgcv.gam = TRUE}.
#' @param permadj a logical variable, whether to adjust the permutated data to cover 
#' the entire range, esp. in case of sparsity. Defaults to \code{FALSE}.
#' @param mc.cores passed to \code{mc.cores} inside of \code{\link[parallel]{mclapply}} 
#' (not available on Windows unless \code{mc.cores = 1}).
#' @param seeds set the seed for the permutation via \code{set.seed(seeds)}
#' @param stepdown a logical variable, which denotes to apply the step-down min P. It defaults to \code{FALSE}.
#' @return 
#' \item{.index}{a vector, evaluation grids.}
#' \item{pval}{a vector or matrix of (adjusted) p values.}
#' \item{vdparam}{a list of paramters of varying distributions.}
#' @note 
#' \itemize{
#' \item If \code{ydata1} and \code{ydata2} are lists of 
#' \code{data.frame}s, the lenghs of two lists must be the same.
#' 
#' \item If \code{mgcv.gam} is \code{TRUE}, \code{...} and \code{exclue} are \code{NULL} 
#' (default settings), then they both default to \cr
#' \code{formula <- list(.value ~ s(.index) + s(.obs, bs = "re"), ~ s(.index))}, \cr
#' \code{family = gauss()} and \code{exclude <- "s(.obs)"}, repectively.
#' 
#' \item Normal distribution (\code{gauss()}) in \code{mgcv::gam} is supported. And \code{\link[gamlss.dist]{gamlss.family}} is supported as well by 
#' \code{DVDtest} for fitting a GAMLSS-type varying distributions with various types of random effect.
#' Note that the permuted data may not match some specific distributions during the permutation. 
#' 
#' }
#' @author Meng Xu, Philip Reiss
#' @references 
#'
#' reiss-EMR18.pdf
#' 
#' Wood, S. N., Pya, N., & Safken, B. (2016). Smoothing parameter and model selection for general 
#' smooth models. Journal of the American Statistical Association, 111(516), 1548-1563.
#' 
#' Rigby, R. A., & Stasinopoulos, D. M. (2005). Generalized additive models for location, scale 
#' and shape. Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(3), 507-554.
#' 
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
#'  
#'  # library(ggplot2)
#'  # ggplot() + geom_line(data = dg1, aes(x = .index,y = .value, group = .obs))
#'  # + geom_line(data = dg2, aes(x = .index, y = .value, group = .obs))
#' 
#'  ngrid <- 50
#'  ev.grid <- seq(0, 1, , ngrid)
#'  nperm. <- 30
#'  
#' ####Estimated with mgcv::gam
#'  library(mgcv)
#'  simu.test1 <- DVDtest(dg1, dg2, nperm. ,ev.grid)
#'  
#' ####Estimated with gamlss::gamlss -! under debugging
#'  \dontrun{
#'  library(gamlss)
#'  library(gamlss.add)
#'  simu.test2 <- DVDtest(dg1, dg2, nperm.,ev.grid, formula = .value ~ ga(~s(.index)) + 
#'                        ga(~s(.obs, bs = "re")), 
#'                        sigma.formula = ~ga(~s(.index)), family = NO, mgcv.gam=FALSE)
#'  }
#' ####Plot
#'  simu.figs <- DVDplot(simu.test1)
#'  simu.figs$pfig
#'  simu.figs$kfig[[1]]
#' 
#' ####Non-normal case -! under debugging
#'  \dontrun{
#'  mu1 <- function(t) 0.2*(p-1)*sin(pi*t)+t+1
#'  mu2 <- function(t) -0.2*(p-1)*sin(pi*t)+t+1
#'  sig1 <- function(t) t+1
#'  sig2 <- sig1
#'  nu1 <- function(t) t+1
#'  nu2 <- nu1
#'  nperson <- 10
#'  library(gamlss)
#'  fun1 <- function(t) rGG(nperson, mu1(t), sig1(t),nu1(t))
#'  fun2 <- function(t) rGG(nperson, mu2(t), sig2(t),nu2(t))
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
#'  
#'  ngrid <- 50
#'  ev.grid <- seq(0, 1, , ngrid)
#'  nperm. <- 30
#'  simu.test3 <- DVDtest(dg1, dg2, nperm.,ev.grid, formula = .value ~ pb(.index), 
#'                        sigma.formula = ~pb(.index),
#'                        nu.formula= ~pb(.index), seeds=123,
#'                        random = ~1|.obs, family = GG, mgcv.gam = FALSE)
#'  }



DVDtest <- function(ydata1, ydata2, nperm, grid, dist.method = "wass", 
                    mgcv.gam = TRUE, ..., exclude = NULL, permadj = FALSE, mc.cores = 1, 
                    seeds = NULL, stepdown = TRUE) {
  argmt <- list(...)
  if (!identical(sapply(ydata1, is.list), sapply(ydata2, is.list))) {
    stop("The types of ydata1 and ydata2 are not same")
  }
  if (!is.list(ydata1[[1]])) {
    ydata1 <- list(ydata1)
    ydata2 <- list(ydata2)
  }
  if (!all(c(c(".obs", ".index", ".value") %in% sapply(ydata1, names), 
             c(".obs", ".index", ".value") %in% sapply(ydata2, names)))) {
    stop("The names of ydata1&2 need to include .index, .obs and .value at least")
  }
  if (length(ydata1) != length(ydata2)) 
    stop("The lengths of ydata1 and ydata2 are not same")
  
  if (mgcv.gam & is.null(exclude) & !is.null(argmt[["formula"]])) 
    stop("The arguments related mgcv::gam are missing, such as formula and exclude")
    
    nroi <- length(ydata1)
    for (k in 1:nroi) {
      ydata1[[k]]$.obs <- factor(ydata1[[k]]$.obs)
      ydata2[[k]]$.obs <- factor(ydata2[[k]]$.obs)
    }
    cat("####perm allocating####", "\n")
    permat.all <- make.perms(ydata1[[1]], ydata2[[1]], nperm, grid, adj = permadj, 
                             seeds = seeds)
    
    cat("####real dist. calculating####", "\n")
    reald <- get_realdist(ydata1 = ydata1, ydata2 = ydata2, 
                          grid = grid, ..., exclude = exclude, dist.method = dist.method, 
                          mc.cores = mc.cores, mgcv.gam = mgcv.gam)
    realdists <- reald$realdists
    vdparam <- reald$vdparam
    
    cat("####perm dist. calculating####", "\n")
    permarray <- wass_perm(nperm = nperm, ydata1 = ydata1, 
                           ydata2 = ydata2, ..., grid = grid, permat = permat.all, exclude = exclude, 
                           mc.cores = mc.cores, dist.method = dist.method, mgcv.gam = mgcv.gam)
    cat("####p-values calculating####", "\n")
    param.array <- get_params(nroi = nroi, nperm = nperm, permarray = permarray, 
                              grid = grid, mc.cores = mc.cores)
    
    # Raw p-values (p.real, p.perm)
    if (!stepdown) {
      p.mat <- get.pval(permarray = permarray, param.array = param.array, 
                        realdists = realdists, nroi = nroi, grid = grid, nperm = nperm)
    } else {
      p.mat <- get.pval.sd(permarray = permarray, param.array = param.array, 
                           realdists = realdists, nroi = nroi, grid = grid, nperm = nperm)
    }
    return(list(.index = grid, pval = p.mat, vdparam = vdparam))
}
