#' @name quantfun
#' @aliases qq
#' @aliases qfun
#' @aliases wqfun
#' 
#' @title Functions to compute (weighted) quantiles
#'
#' @description
#' This set of functions allows you to compute weighted quantiles.
#' \code{qq} returns the quantile from a sequence of values with associated
#' probabilities. \code{qfun} computes unweighted quantiles using \code{qq},
#' and \code{wqfun} computes the weighted quantiles given some weight for
#' the different values. \code{quantfun} computes weighted quantiles for an
#' array of input data with equal weight to rows (and allowing for missing
#' values in the columns).
#'
#' @param x numerical input (vector or matrix for \code{quantfun})
#' @param initw initial weight for individual rows in \code{x}
#' @param quants multiple quantiles to be computed
#' @param type Assumptions on how quantiles are calculated (see details)
#' 
#' @details
#' The different \code{type}s govern the probabilities associated with 
#' the sequence of ordered values. \code{type = 1} sets the probability of 
#' the minimum to \code{0.5/n} where \code{n} is the number of values in 
#' \code{x}, whereas \code{type = 2} sets the probability of the minimum
#' to \code{0}. In both cases the 0th and 100th percentile correspond to the
#' minimum and maximum of the values in \code{x}.
#' 
#' @keywords utilities
#' 
#' @rdname quantfun
#' @export
quantfun <- function(x, initw=1, quants=c(0.1, 0.9, 0.5), type=1){
  initw <- rep(initw, length=length(x))
  if(sum(!is.na(x)) < 2) {
    xout <- rep(NA, length(quants))
  } else if (length(dim(x)) == 0){
    xout <- wqfun(x, weights=initw, p=quants, type=type)
  } else {
    xw <- initw * rep(1/apply(!is.na(as.matrix(x)), 1, sum), length=length(x))
    xout <- wqfun(x, weights=xw, p=quants)
    
  }
  return(xout)
}
#'
#' @rdname quantfun
#' @param px probability of values in \code{x}
#' @param p quantile to be computed computed
#' @export

qq <- function(x, px, p){
  k <- as.numeric(cut(p, px, include.lowest=TRUE))
  out <- x[k] + (p - px[k]) / (px[pmin(k + 1, length(x))] - px[k]) * (x[pmin(k+1, length(x))] - x[k])
  out[p > max(px)] <- max(x)
  out[p < min(px)] <- min(x)
  return(out)
}  

#' @rdname quantfun
#' @export
qfun <- function(x, p, type=1){
  x <- sort(x[!is.na(x)])
  if (type == 1){
    px <-  (seq(along=x) - 0.5)/length(x)
  } else {
    px <- seq(0,1,length=length(x))
  }  
  out <- qq(x, px, p)
  return(out)
}

#' @rdname quantfun
#' @param weights weights for the values in \code{x}
#' @export
wqfun <- function(x, weights, p, type=1){
  if (sum(!is.na(x)) >= 2){
    ix <- sort(x, ind=T)$ix
    xx <- x[!is.na(x)][ix]
    xw <- weights[!is.na(x)][ix]
    xw <- xw / sum(xw)
    ## remove all entries with zero weight
    xx <- xx[xw > .Machine$double.neg.eps]
    xw <- xw[xw > .Machine$double.neg.eps]
    px <- 1/sum(xw) * (cumsum(xw) - xw/2)
    if (type == 2){
      px <- (px - px[1]) / (px[length(px)] - px[1])
    }   
    out <- qq(xx, px, p)
  } else {
    out <- rep(NA, length(p))
  }
  return(out)
}

