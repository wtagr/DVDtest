#' Difference between Varying Distributions Test (DVDtest)
#' 
#' Testing the difference of two varying distributions.
#' 
#' 
#' 
#' @param ydata1 a \code{data.frame} or a \code{list} of \code{data.frame}, containing 
#' at least 3 columns called \code{.obs}, \code{.index} and \code{.value} which 
#' specify which curve the point belongs to (\code{.obs}) at which (\code{.index}) 
#' it was observed and the observed value (\code{.value}). See details in the package 
#' \code{refund}. Other columns are available as well for modeling the varying 
#' distributions.
#' @param ydata2 same as \code{ydata1}. 
#' @param nperm a scalar, number of permutation
#' @param grid a vector, evaluation grids of \code{.index}
#' @param dist.method the distance measure to be used. This must be one of Wasserstein 
#' distance (\code{'wass'}), \code{'L2'} distance, \code{'L1'} distance and \code{'Hellinger'}. 
#' It defaults to \code{'wass'}.
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
#' @param seeds a scalar, to set the seed for the permutation via \code{set.seed(seeds)}. It defaults to \code{NULL}.
#' @param chunksize a scalar, number of permutations per core/threat when using \code{mcapply}. It defaults to 2.
#' @param report.every a scalar, to report the current number of permutations.
#' @param stepdown a logical variable, which denotes to apply the step-down min P.
#' @return 
#' \item{.index}{a vector, evaluation grids.}
#' \item{pval}{a vector or matrix of (adjusted) p values.}
#' \item{vdparam}{a list of paramters of varying distributions.}
#' @note 
#' \itemize{
#' \item If \code{ydata1} and \code{ydata2} are lists of 
#' \code{data.frame}s, the lengths of two lists must be the same.
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
#' @import mgcv gamlss gamlss.dist parallel
#' @examples
#' 
#' library(DVDtest)
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
#'  # library(ggplot2)
#'  # ggplot() + geom_line(data = dg1, aes(x = .index,y = .value, group = .obs))
#'  # + geom_line(data = dg2, aes(x = .index, y = .value, group = .obs))
#' 
#'  ngrid <- 50
#'  ev.grid <- seq(0, 1, , ngrid)
#'  nperm. <- 40
#'  
#' ####Estimated with mgcv::gam
#'  library(mgcv)
#'  simu.test1 <- DVDtest(dg1, dg2, nperm., ev.grid)
#'  
#' ####Estimated with gamlss::gamlss
#'  library(gamlss)
#'  simu.test2 <- DVDtest(dg1, dg2, nperm.,ev.grid, formula = .value ~ pb(.index), 
#'                sigma.formula = ~pb(.index), random = ~1|.obs, family = NO, mgcv.gam=FALSE)
#'                
#' ####Plot
#'  simu.figs <- DVDplot(simu.test1)
#'  simu.figs$pfig
#'  simu.figs$kfig[[1]]
#'  
#' 
#' ####Non-normal case
#' \dontrun{
#'  p <- 6
#'  mu1 <- function(t) 0.2*(p-1)*sin(pi*t)+t+1
#'  mu2 <- function(t) -0.2*(p-1)*sin(pi*t)+t+1
#'  sig1 <- function(t) t+1
#'  sig2 <- sig1
#'  nu1 <- function(t) t+1
#'  nu2 <- nu1
#'  nperson <- 10
#'  
#'  fun1 <- function(t) rGG(nperson, mu1(t), sig1(t),nu1(t))
#'  fun2 <- function(t) rGG(nperson, mu2(t), sig2(t),nu2(t))
#'  tp <- seq(0,1,,10)
#'  data1 <- sapply(tp,fun1)
#'  data2 <- sapply(tp,fun2)
#'  
#'  library(reshape2)
#'  colnames(data2) <- colnames(data1) <- tp
#'  rownames(data2) <- 1:nperson+nperson
#'  dg1 <- melt(data1)
#'  dg2 <- melt(data2)
#'  colnames(dg1) <- colnames(dg2) <- c('.obs','.index','.value')
#'  
#'  ngrid <- 50
#'  ev.grid <- seq(0, 1, , ngrid)
#'  nperm. <- 40
#'  
#'  simu.test3 <- DVDtest(dg1, dg2, nperm.,ev.grid, formula = .value ~ pb(.index), 
#'                        sigma.formula = ~pb(.index),
#'                        nu.formula= ~pb(.index), seeds=123,
#'                        random = ~1|.obs, family = GG, mgcv.gam = FALSE)
#'}



DVDtest <- function(ydata1, ydata2, nperm = 60, grid, dist.method = "wass", 
                    mgcv.gam = TRUE, ..., exclude = NULL, permadj = FALSE, mc.cores = 1, 
                    chunksize = 2, seeds = NULL, report.every = 10, stepdown = FALSE) {
  argmt <- list(...)
  if (!identical(sapply(ydata1, is.list), sapply(ydata2, is.list))) {
    stop("The types of ydata1 and ydata2 are not same")
  }
  if (mc.cores < 1L) stop("'mc.cores' must be >= 1")
  os <- Sys.info()
  if ("Windows"%in%os ) {
    if (mc.cores > 1) warning( "'mc.cores' > 1 is not supported on Windows, and it force to set as 'mc.cores = 1'" )
    mc.cores <- 1
  }

  if (nperm/chunksize/mc.cores-as.integer(nperm/chunksize/mc.cores) > 1e-12)
    stop("We suggest the arguments 'chunksize' and 'mc.cores' are multiply factors of number of permutation 'nperm'.")
  
  if (!is.list(ydata1[[1]])) {
    ydata1 <- list(ydata1)
    ydata2 <- list(ydata2)
  } 
  if (!all(c(c(".obs", ".index", ".value") %in% sapply(ydata1, names), 
             c(".obs", ".index", ".value") %in% sapply(ydata2, names)))) {
    stop("The names of ydata1&2 need to include .index, .obs and .value")
  }
  if (length(ydata1) != length(ydata2)) 
    stop("The length of ydata1 and ydata2 are not same")
  
  if (mgcv.gam & is.null(exclude) & !is.null(argmt[["formula"]])) 
    stop("The arguments related mgcv::gam are missing, such as formula and exclude")
  
  nroi <- length(ydata1)
  for (k in 1:nroi) {
    ydata1[[k]]$.obs <- factor(ydata1[[k]]$.obs)
    ydata2[[k]]$.obs <- factor(ydata2[[k]]$.obs)
  }
  
  permat.all <- make.perms(ydata1[[1]], ydata2[[1]], nperm, grid, adj = permadj, 
                           seeds = seeds)
  cat("#### perm allocated ####", "\n")
  
  reald <- get_realdist(ydata1 = ydata1, ydata2 = ydata2, grid = grid, 
                        ..., exclude = exclude, dist.method = dist.method, mc.cores = mc.cores, mgcv.gam = mgcv.gam)
  realdists <- reald$realdists
  cat("#### real dist. calculated ####", "\n")
  
  vdparam <- reald$vdparam
  
  permatlist <- list()
  # number of permutations per core
  for (j in 1:(nperm/chunksize)) {
    mm <- chunksize * (j - 1) + 1
    MM <- chunksize * j
    permatlist[[j]] <- permat.all[mm:MM, ]
  }

  perm.chunk <- function(whichmat) {
    ob <- list()
    for (k in 1:nroi) {
      ob[[k]] <- wass.perm(ydata1[[k]], ydata2[[k]], permat = permatlist[[whichmat]], 
                           grid = grid, dist.method = dist.method, mgcv.gam = mgcv.gam)
    }
    ob
  }

  ptemp = list()
  for (j in 1:(nperm/chunksize/mc.cores)) {
    if (!j%%report.every) cat("*** Perm", j*chunksize*mc.cores, "***\n")
    lb <- mc.cores * (j - 1) + 1
    ub <- mc.cores * j
    ptemp[[j]] = mclapply(lb:ub, perm.chunk, mc.cores = mc.cores, mc.preschedule = FALSE)
    # save(list=paste0('p',j), file=paste0('perms/p',j,'.Rdata'))
  }
  cat("*** Perm", j*chunksize*mc.cores, "***\n")
  # Make a 101 x 74 x nperm. array of permdists nerrvec <- rep(0,nroi)
  ndone <- nperm/chunksize/mc.cores

  permarray <- array(dim = c(nperm, length(grid), nroi))
  pcount <- 0
  for (k in 1:ndone) {
    pp <- ptemp[[k]]
    for (co in 1:mc.cores) {
      lob <- pcount + 1
      pcount <- pcount + chunksize
      for (reg in 1:nroi) {
        permarray[lob:pcount, , reg] <- pp[[co]][[reg]]$permdists
        # nerrvec[reg] <- nerrvec[reg] + pp[[co]][[reg]]$nerror
      }
    }
  }
  
  # permarray <- wass_perm(vdFun = vdFun, nperm = nperm, ydata1 = ydata1,
  # ydata2 = ydata2,..., grid = grid, permat = permat.all, exclude =
  # exclude, mc.cores = mc.cores, dist.method = dist.method)
  cat("#### perm dist. calculated ####", "\n")
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
