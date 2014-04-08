\name{open.ncdf}
\alias{open.ncdf}
\title{New \code{open.ncdf} for CF compatibility}
\usage{
\method{open}{ncdf}(con, write = FALSE, readunlim = TRUE, verbose = FALSE,
  ...)
}
\arguments{
  \item{con}{Name of the existing netCDF file to be
  opened.}

  \item{write}{If \code{FALSE} (default), then the file is
  opened read-only. If \code{TRUE}, then writing to the
  file is allowed.}

  \item{readunlim}{When invoked, this function reads in the
  values of all dimensions from the associated variables.
  This can be slow for a large file with a long unlimited
  dimension. If set to \code{FALSE}, the values for the
  unlimited dimension are not automatically read in (they
  can be read in later, manually, using
  \code{\link[ncdf]{get.var.ncdf}()}).}

  \item{...}{Arguments passed to or from other methods.}

  \item{verbose}{If \code{TRUE}, then messages are printed
  out during execution of this function.}
}
\description{
This function replaces the \code{\link[ncdf]{open.ncdf}} in
the \code{\link{ncdf}} package to ensure compliance with CF
standard with regard to the missing value attribute.
}

