% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.params.R
\name{get.params}
\alias{get.params}
\title{Fitting the gamma distribution (Pearson type III) for each k}
\usage{
get.params(k, nperm, permarray, grid)
}
\arguments{
\item{k}{a scalar, \code{k}th data.frame in \code{ydata1&2}}

\item{nperm}{a scalar, number of permutation}

\item{permarray}{an array, permuted-data distances from \code{wass_perm}}

\item{grid}{a vector, evaluation grids of \code{.index}}
}
\value{
a list of mean, sigma, nu and the AIC
}
\description{
Fitting by the permuted-data distances
}
\seealso{
\code{\link{get_params}}
}
\author{
Philip Reiss, Meng Xu
}
\keyword{internal}
