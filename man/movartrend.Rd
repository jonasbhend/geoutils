\name{movartrend}
\alias{artrend}
\alias{movartrend}
\alias{trendcov}
\alias{trendresid}
\title{Compute trends on NetCDF objects}
\usage{
movartrend(data, trndln, x = NULL, missthresh = 0.8, robust = FALSE,
  gaps = FALSE, plotable = TRUE, resid = FALSE)

trendcov(data, x = NULL, missthresh = 0.8)

trendresid(data, x = NULL, missthresh = 0.8)

artrend(data, x = NULL, missthresh = 0.8)
}
\arguments{
  \item{data}{object of type NetCDF containing the time
  series}

  \item{trndln}{length of trend to be computed (in years if
  time attribute is present)}

  \item{x}{optional covariate to be used as predictor}

  \item{missthresh}{threshold to indicate what fraction of
  non-missing values is at least required to run the
  regression}

  \item{robust}{logical, should robust regression be used?}

  \item{gaps}{logical, should gaps be allowed for?}

  \item{plotable}{logical, should output be rearrange to
  make it plotable?}

  \item{resid}{logical, should residuals of regression be
  returned?}
}
\description{
Compute trends on NetCDF objects where \code{movartrend}
computes moving linear regressions on the time series,
\code{trendresid} and \code{trendcov} compute simple linear
regression with the latter returning the covariance matrix
of the standard errors (i.e. for multivariate analyses).
}
\keyword{utilities}

