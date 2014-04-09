#' @name colourramp
#' @aliases .getval .colseq
#' 
#' @title Colour palettes
#' 
#' @description 
#' Functions to help set up colour palettes for use in plots.
#' 
#' @param n number of colour values to be returned
#' @param ramp text specifying what ramp to use
#' @param start darkness for darkest colour (0 is black)
#' @param log logical, should interval be divided linearly?
#' 
#' @keywords utilities
#' @export
colourramp <- function(n, ramp='redblue', start=0.4, log=F){
  numx <- floor(n/2)
  uppercol <- if (ramp == 'redblue') .redwhite else if (ramp == 'bluebrown') .bluewhite else if (ramp == 'greenbrown') .greenwhite else stop('Colour ramp not implemented')
  undercol <- if(ramp == 'redblue') .bluewhite else if (ramp == 'bluebrown') .brownwhite else if (ramp == 'greenbrown') .brownwhite else stop('Colour ramp not implemented')
  upperseq <- .colseq((numx + 1), uppercol, start=start, log=log)
  underseq <- .colseq((numx + 1), undercol, start=start, log=log)
  if (2*numx == n){
    colour <- c(underseq[1:numx], rev(upperseq[1:numx]))
  } else {
    colour <- c(underseq, rev(upperseq[1:numx])) 
  }
  return(colour)
}


#' @rdname colourramp
#' @param ind An index specifying the spacing of the values (0 - 1)
#' @param val The values to which a smoothed spline is fitted
#' @param predind The location for which the smoothed spline function
#'           is evaluated (0 - 1)
#' @param smooth the smoothing as \code{spar} in \code{\link{smooth.spline}}
#' 
#' @export
.getval <- function(ind, val, predind, smooth) {
  # Function to get smoothed values in the interval 0 - 1    
  fitfun  <- smooth.spline(ind, val, spar=smooth)
  fitval  <- predict(fitfun, predind)$y
  if (min(fitval) < 0) fitval <- fitval - min(fitval)
  if (max(fitval) > 1) fitval <- fitval / max(fitval)
  fitval
}

#' @rdname colourramp
#' @param x Number of colors 
#' @param tmp Array of rgb values (dim=c(ncolors,3)) along which
#'            the color sequence is interpolated
#' @param stop the upper bound for brightness of colours (1 is white)
#' @export
.colseq <- function(x, tmp, start=0, stop=1, log=F, smooth=0.4){
  # .colseq computes a color sequence of x rgb colors along the indicated
  # colors by tmp (Attention: 'ind' is set to .spacing to match the 
  # colours provided by .bluewhite, etc.)  
  ind <- seq(0,1,length=dim(tmp)[1])
  predind <- seq(start, stop, length=x)
  if (log & min(predind) > 0){
    predind <- log(predind)
    mini    <- min(predind)
    maxi    <- max(predind)
    predind <- (predind - mini)/(maxi-mini)*(stop-start) + start
  }
  
  r       <- .getval(ind, tmp[,1], predind, smooth=smooth)
  g       <- .getval(ind, tmp[,2], predind, smooth=smooth)
  b       <- .getval(ind, tmp[,3], predind, smooth=smooth)
  
  colour  <- rgb(r,g,b)
  colour
}

#' Colours used in the geoutils package
#' 
#' Provides a few (hidden) colour sequences for use with colfun
#' 
#' @docType data
#' @keywords datasets
#' @format matrices with RGB colour values (hex colours for \code{.soil}
#'         and \code{.water})
#' @aliases .hsurf .gpcc .redwhite .bluewhite .brownwhite .greenwhite
#' @name geocolours
NULL
