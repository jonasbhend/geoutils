\name{plot.NetCDF}
\alias{plot.NetCDF}
\title{Function to plot objects of class 'NetCDF'}
\usage{
\method{plot}{NetCDF}(x, type = if (nrow(x) == 1) "ts" else "mean", ti = if
  ("time" \%in\% names(attributes(x))) seq(along = attr(x, "time")) else 1,
  si = 1, levs = NULL, col = NULL, pt.cex = 2, symmetric = FALSE,
  cut = TRUE, add = FALSE, xlab = "", ylab = "", xlim = NULL,
  ylim = NULL, xaxs = "i", yaxs = "i", lty = seq(along = si), lwd = 3,
  nlevs = 12, colramp = "redblue", seas = F, ...)
}
\arguments{
  \item{x}{input object of class NetCDF (max. 2-d)}

  \item{ti}{index of time step to be plotted (map only)}

  \item{si}{index of spatial grid to be plotted}

  \item{levs}{levels to be used for contouring}

  \item{col}{colours to be used for contouring or lines}

  \item{pt.cex}{expansion factor for points}

  \item{type}{of plot ('ts' for lines, 'mean' and 'trend'
  for maps)}

  \item{symmetric}{logical, should contouring be symmetric
  about zero}

  \item{cut}{logical, should contouring be cut to range of
  data}

  \item{add}{logical, should plot be added}

  \item{xlab,ylab}{axes labels}

  \item{xlim,ylim}{range of axes}

  \item{xaxs,yaxs}{type of axes}

  \item{lty}{line type for line plots}

  \item{lwd}{line width for line plots}

  \item{nlevs}{number of contour interval (input to
  \code{pretty})}

  \item{colramp}{colour ramp to be used with colourramp
  function (default)}

  \item{seas}{logical, should seasons be plotted
  individually?}
}
\description{
This function allows plotting line plots and maps of NetCDF
objects using the respective attributes.
}
\examples{
x <- array(rnorm(10000), c(100,100)) + as.vector(outer(sort(rnorm(10)), sort(rnorm(10, mean=15)), '+')) + rep(c(-2, 0, 3, 1), each=100)
mostattributes(x) <- list(class='NetCDF', time=seq(1950.08333, 1975, 0.25), dim=c(100,100), lon=rep(1:10, 10), lat=rep(1:10, each=10))
plot(x, type='ts', si=1:3)

## also plot a map (very rudimentary)
plot(x, type='trend', pt.cex=3)
}
\keyword{plot}

