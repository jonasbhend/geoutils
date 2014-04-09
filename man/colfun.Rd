\name{colfun}
\alias{.colseq}
\alias{.getval}
\alias{bbfun}
\alias{colfun}
\alias{gbfun}
\alias{rbfun}
\title{Colour palettes}
\usage{
.getval(ind, val, predind, smooth)

.colseq(x, tmp, start = 0, stop = 1, log = F, smooth = 0.4)

rbfun(n, start = 0.4, log = F)

gbfun(n, start = 0.4, log = F)

bbfun(n, start = 0.4, log = F)
}
\arguments{
  \item{ind}{An index specifying the spacing of the values
  (0 - 1)}

  \item{val}{The values to which a smoothed spline is
  fitted}

  \item{predind}{The location for which the smoothed spline
  function is evaluated (0 - 1)}

  \item{smooth}{the smoothing as \code{spar} in
  \code{\link{smooth.spline}}}

  \item{x}{Number of colors}

  \item{tmp}{Array of rgb values (dim=c(ncolors,3)) along
  which the color sequence is interpolated}

  \item{start,stop}{the range for which colours should be
  interpolated (0 - 1)}

  \item{log}{logical, should interval between \code{start}
  and \code{stop} be divided linearly?}

  \item{n}{number of colours to return}
}
\description{
Functions to help set up colour palettes for use in plots.
}
\keyword{utilities}

