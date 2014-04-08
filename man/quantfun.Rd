\name{quantfun}
\alias{qfun}
\alias{qq}
\alias{quantfun}
\alias{wqfun}
\title{Functions to compute (weighted) quantiles}
\usage{
quantfun(x, initw = 1, quants = c(0.1, 0.9, 0.5), type = 1)

qq(x, px, p)

qfun(x, p, type = 1)

wqfun(x, weights, p, type = 1)
}
\arguments{
  \item{x}{numerical input (vector or matrix for
  \code{quantfun})}

  \item{initw}{initial weight for individual rows in
  \code{x}}

  \item{quants}{multiple quantiles to be computed}

  \item{type}{Assumptions on how quantiles are calculated
  (see details)}

  \item{px}{probability of values in \code{x}}

  \item{p}{quantile to be computed computed}

  \item{weights}{weights for the values in \code{x}}
}
\description{
This set of functions allows you to compute weighted
quantiles. \code{qq} returns the quantile from a sequence
of values with associated probabilities. \code{qfun}
computes unweighted quantiles using \code{qq}, and
\code{wqfun} computes the weighted quantiles given some
weight for the different values. \code{quantfun} computes
weighted quantiles for an array of input data with equal
weight to rows (and allowing for missing values in the
columns).
}
\details{
The different \code{type}s govern the probabilities
associated with the sequence of ordered values. \code{type
= 1} sets the probability of the minimum to \code{0.5/n}
where \code{n} is the number of values in \code{x}, whereas
\code{type = 2} sets the probability of the minimum to
\code{0}. In both cases the 0th and 100th percentile
correspond to the minimum and maximum of the values in
\code{x}.
}
\keyword{utilities}

