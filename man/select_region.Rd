\name{select_region}
\alias{select_region}
\title{Select a region}
\usage{
select_region(x, regi = 1)
}
\arguments{
  \item{x}{input array of class NetCDF}

  \item{regi}{index or name of region to select}
}
\description{
Function to select regions (2nd last dimension in NetCDF
array).
}
\details{
Named regions to select only work if the NetCDF object
contains an attribute 'regions' that contains the region
names referred to in \code{regi}.
}
\keyword{utilities}

