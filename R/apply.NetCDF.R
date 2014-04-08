#' Apply a function to specified dimensions of an 'NetCDF' object
#'
#' This function works as \code{\link{apply}} but also applies functions
#' only to time steps with the same season/month, so that seasonal analysis
#' can be easily performed.
#' 
#' @param x input array of class 'NetCDF'
#' @param dims dimensions that are retained in the output (see \code{\link{apply}})
#' @param fun function to be applied
#' @param ... additional arguments to \code{fun}
#' @details 
#' Although apply.NetCDF suggest to be a generic function that can be
#' applied on objects of class 'NetCDF' with \code{apply(x, ..)}, this
#' is deliberately not the case. \code{\link{apply}} is not a generic function
#' and therefore you will have to call \code{apply.NetCDF} specifically
#' to call this function.
#' @keywords utilities
#' @export
apply.NetCDF <- function(x, dims, fun, ...){
  ## applies a function to a .NetCDF object
  ## taking into account a time attribute if present and rearranging
  ## the array accordingly
  if (is.null(attr(x, 'time')) | length(dim(x)) %in% dims){
    xout <- apply(x, dims, fun, ...)
    for (atn in setdiff(names(attributes(x)), c('dim', 'dimnames'))){
      attr(xout, atn) <- attr(x, atn)
    }
  } else {
    ## deal with the case of applying function to time dimension
    nseas <- median(table(ceiling(attr(x, 'time'))))
    indims <- dim(x)
    tmpdims <- c(dim(x)[-length(dim(x))], nseas, dim(x)[length(dim(x))]/nseas)
    outdims <- c(indims[dims], nseas)
    ndims <- c(dims, length(dim(x)))
    xtmp <- array(x, tmpdims)
    xout <- apply(xtmp, ndims, fun, ...)
    ## check whether output has same length as input (reshuffle if yes)
    if (length(xout) == length(x)){
      xout <- array(aperm(xout, c(2:length(dim(xout)), 1)), dim(x))
    }
    ## rewrite attributes
    for (atn in setdiff(names(attributes(x)), c('dim', 'dimnames'))){
      attr(xout, atn) <- attr(x, atn)
    }
    if (length(xout) != length(x)){
      ## redefine time attribute
      attr(xout, 'time') <- attr(x, 'time')[1:nseas]
      attr(xout, 'time_average') <- range(floor(attr(x, 'time')))
    }
  }
  return(xout)
}
