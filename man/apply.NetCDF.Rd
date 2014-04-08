\name{apply.NetCDF}
\alias{apply.NetCDF}
\title{Apply a function to specified dimensions of an 'NetCDF' object}
\usage{
apply.NetCDF(x, dims, fun, ...)
}
\arguments{
  \item{x}{input array of class 'NetCDF'}

  \item{dims}{dimensions that are retained in the output
  (see \code{\link{apply}})}

  \item{fun}{function to be applied}

  \item{...}{additional arguments to \code{fun}}
}
\description{
This function works as \code{\link{apply}} but also applies
functions only to time steps with the same season/month, so
that seasonal analysis can be easily performed.
}
\details{
Although apply.NetCDF suggest to be a generic function that
can be applied on objects of class 'NetCDF' with
\code{apply(x, ..)}, this is deliberately not the case.
\code{\link{apply}} is not a generic function and therefore
you will have to call \code{apply.NetCDF} specifically to
call this function.
}
\keyword{utilities}

