\name{detect_break}
\alias{detect_break}
\title{Detect step changes}
\usage{
detect_break(x, breaktime = 2005, n = 50, detrend = TRUE)
}
\arguments{
  \item{x}{input object of type NetCDF}

  \item{breaktime}{last time step (year) of first period}

  \item{n}{number of years in climatology (before and after
  break point)}

  \item{detrend}{logical, should linear trend be accounted
  for?}
}
\description{
This function checks whether there are step changes present
in a time series at a pre-specified time (generally 2005
for CMIP5 historical to RCP scenario transition).
}
\keyword{utilities}

