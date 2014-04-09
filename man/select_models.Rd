\name{select_models}
\alias{select_ensembles}
\alias{select_models}
\title{Select models and ensemble members}
\usage{
select_models(x, modi = 1:nrow(x))

select_ensembles(x, ensi = 1:ncol(x))
}
\arguments{
  \item{x}{object of class 'NetCDF'}

  \item{modi}{indices of models to select (first dimension
  in x)}

  \item{ensi}{indices (or text) of runs to select (second
  dimension of x)}
}
\description{
Functions to select models \code{select_models} and
ensemble members \code{select_ensembles} from an array of
model results.
}
\keyword{utilities}

