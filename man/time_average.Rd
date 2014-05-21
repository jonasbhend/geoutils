\name{time_average}
\alias{time_average}
\title{Create n-year averages of an NetCDF object}
\usage{
time_average(x, n, offset = NULL, type = "start")
}
\arguments{
  \item{x}{array of class 'NetCDF'}

  \item{n}{number of time steps to average over}

  \item{offset}{offset from which to start (in years)}

  \item{type}{Alternative to specifying offset, one of
  start (the default), even or end, see details}
}
\description{
This function computes n-year averages from an existing
NetCDF time series read in with \code{\link{readNetCDF}}
}
\details{
The automated way to set the offset works as following.
\code{start} starts the \code{n}-year averages at the
beginning of the time series, \code{end} sets the offset so
that the last average ends in the last year of the time
series, \code{even} sets the offset so that
\code{modulo(year, n) == 0} for the first year in the first
average (i.e. for \code{n=5} averages such as 1910-1914)
}
\keyword{utilities}

