% Generated by roxygen2 (4.0.1): do not edit by hand
\name{vert.poly}
\alias{vert.poly}
\title{Plot shaded area between lines}
\usage{
vert.poly(y1, y2, x = seq(along = y1), horizontal = TRUE, ...)
}
\arguments{
\item{y1,y2}{vectors for vertical lines}

\item{x}{vector with x intersects (optional)}

\item{horizontal}{logical, x-y axis inverted if \code{FALSE}}

\item{...}{additional arguments passed to \code{\link{polygon}}}
}
\description{
Function to simplify plotting shaded areas between two lines.
Functionality to deal with missing values is added (breaks in
the shaded areas where there are missing values).
}
\keyword{plot}

