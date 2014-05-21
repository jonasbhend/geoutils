\name{anomalies}
\alias{anomalies}
\title{Compute anomalies with respect to reference period}
\usage{
anomalies(x, anomaly = c(1961, 1990), relative = FALSE, na.thresh = 0.8,
  rel.thresh = 1e+06, zero.thresh = 1)
}
\arguments{
  \item{x}{object of type NetCDF}

  \item{anomaly}{start and end year of reference period}

  \item{relative}{logical, should relative anomalies (in
  percent) be computed?}

  \item{na.thresh}{minimum fraction of non-missing values
  in reference period}

  \item{rel.thresh}{maximum relative anomalies that get
  computed}

  \item{zero.thresh}{maximum number of zeros in reference
  for relative anomalies to be computed}
}
\description{
This function computes anomalies with respect to a given
reference period. For the function to work, the reference
period must be fully in the time series provided.
}
\keyword{utilities}

