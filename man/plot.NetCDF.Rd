% Generated by roxygen2 (4.0.1): do not edit by hand
\name{plot.NetCDF}
\alias{plot.NetCDF}
\title{Function to plot objects of class 'NetCDF'}
\usage{
\method{plot}{NetCDF}(x, type = if (nrow(x) == 1) "ts" else "mean",
  ti = 1:dim(x)[length(dim(x))], si = 1, levs = NULL, col = NULL,
  na.col = NULL, pt.cex = NULL, symmetric = FALSE, cut = TRUE,
  add = FALSE, xlab = "", ylab = "", xlim = NULL, ylim = NULL,
  xaxs = "i", yaxs = "i", lty = seq(along = si), lwd = 3, nlevs = 12,
  colramp = "redblue", seas = F, ...)
}
\arguments{
\item{x}{input object of class NetCDF (max. 2-d)}

\item{type}{of plot ('ts' for lines, 'mean' and 'trend' for maps)}

\item{ti}{index of time step to be plotted (map only)}

\item{si}{index of spatial grid to be plotted}

\item{levs}{levels to be used for contouring}

\item{col}{colours to be used for contouring or lines}

\item{na.col}{colour to be used for missing values (spatial plot only)}

\item{pt.cex}{expansion factor for points}

\item{symmetric}{logical, should contouring be symmetric about zero}

\item{cut}{logical, should contouring be cut to range of data}

\item{add}{logical, should plot be added}

\item{xlab,ylab}{axes labels}

\item{xlim,ylim}{range of axes}

\item{xaxs,yaxs}{type of axes}

\item{lty}{line type for line plots}

\item{lwd}{line width for line plots}

\item{nlevs}{number of contour interval (input to \code{pretty})}

\item{colramp}{colour ramp to be used with colourramp function (default)}

\item{seas}{logical, should seasons be plotted individually?}

\item{...}{additional arguments passed to \code{plot.default}}
}
\description{
This function allows plotting line plots and maps of NetCDF objects using
the respective attributes.
}
\examples{
tas <- readNetCDF(system.file("extdata", "annual_CRUTEMv3_1961-90.nc", package="geoutils"), varname="temp")
## average over all time steps
plot(tas)

## also plot the trends
plot(tas, type='trend', sym=TRUE)

## plot the same with legend
layout(matrix(1:2, 2, 1), height=c(5, lcm(2.5)))
par(mar=c(2,2,1,1))
tmp <- plot(tas, type='trend', sym=TRUE, na.col='lightgrey')
## if pt.cex is not set, the default pt.cex will be displayed
map(add=T, interior=F)
plot_colourbar(tmp, units='K/a')
}
\keyword{plot}

