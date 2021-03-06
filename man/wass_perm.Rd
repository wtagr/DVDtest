% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wass_perm.R
\name{wass_perm}
\alias{wass_perm}
\title{Calculating the distances via permuted data}
\usage{
wass_perm(
  nperm,
  ydata1,
  ydata2,
  grid,
  permat,
  ...,
  exclude,
  mc.cores,
  dist.method,
  mgcv.gam
)
}
\arguments{
\item{nperm}{a scalar, number of permutation}

\item{ydata1}{see \code{DVDtest}}

\item{ydata2}{see \code{DVDtest}}

\item{grid}{see \code{\link{DVDtest}}}

\item{permat}{a result of \code{make.perm}}

\item{\dots}{partial arguments of \code{vdFun}}

\item{exclude}{an argument of \code{predict}}

\item{mc.cores}{a scalar, an argument of \code{mclapply}}

\item{dist.method}{see \code{DVDtest}}
}
\value{
an array of distances
}
\description{
Calculating the distances via permuted data
}
\seealso{
\code{wass.perm}
}
\author{
Meng Xu, Philip Reiss
}
\keyword{internal}
