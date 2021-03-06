% Generated by roxygen2 (4.0.1): do not edit by hand
\name{fdrfun}
\alias{fdrfun}
\title{Compute p-values when controlling the False Discovery Rate}
\usage{
fdrfun(p, q = 0.05, log = TRUE, tight = TRUE)
}
\arguments{
\item{p}{p-values of multiple tests}

\item{q}{probability for which FDR should be controlled}

\item{log}{logical, are p-values log-transformed?}

\item{tight}{logical, should tighter FDR control be used}
}
\description{
Function to compute adjusted critical p-values and the fraction of
test for which the alternative hypothesis might be true according
to Ventura et al. (2004).
}
\references{
Ventura, V., C.J. Paciorek, and J.S. Risbey (2004). Controlling the proportion of falsely rejected hypotheses when conducting multiple tests with climatological data. Journal of Climate, 17 (22).
}
\keyword{utilities}

