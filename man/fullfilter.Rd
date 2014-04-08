\name{fullfilter}
\alias{applyfilter}
\alias{desired_response}
\alias{fullfilter}
\alias{get_filter_weights}
\alias{get_response}
\alias{make_filter}
\title{Filter functions for 'NetCDF' objects}
\usage{
fullfilter(x, h)

applyfilter(x, h)

get_response(h, accuracy = 1e-05)

desired_response(fr, type = "lowpass", accuracy = 1e-05)

get_filter_weights(n, fr, accuracy = 1e-05)

make_filter(n, fr, type = "lowpass", ...)
}
\arguments{
  \item{x}{time series or 'NetCDF' object}

  \item{h}{filter weights}

  \item{accuracy}{numerical accuracy to which the response
  should be computed}

  \item{fr}{cutoff frequencies}

  \item{type}{one of lowpass, bandpass, or highpass}

  \item{n}{number of time steps used in digital filter}

  \item{...}{additional arguments passed to
  \code{get_filter_weights}}
}
\description{
\code{fullfilter} applies a moving average filter to the
full time series (no missing values, decreasing filter
length at start and end). \code{applyfilter} applies
\code{fullfilter} to the time dimension of an object of
class 'NetCDF'. \code{make_filter} helps you set up a low-
band- or highpass filter with the desired response
specified in \code{desired_response}.
\code{get_filter_weights} provides the weights given a
cutoff frequency. \code{get_response} can be used to
illustrate and analyse the filter response.
}
\examples{
xx <- rnorm(100) + sin(seq(0,4*pi, length=100))
plot(xx, type='l')
hh <- make_filter(31, 10, type='lowpass')
lines(filter(xx, hh), lwd=2)
plot(get_response(hh), type='l')
abline(v=1/10, lty=2, lwd=2)
}
\keyword{utilities}

