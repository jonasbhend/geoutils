#' @name colfun
#' @aliases .getval .colseq rbfun gbfun bbfun
#' 
#' @title Colour palettes
#' 
#' @description 
#' Functions to help set up colour palettes for use in plots.
#' 
#' @param ind An index specifying the spacing of the values (0 - 1)
#' @param val The values to which a smoothed spline is fitted
#' @param predind The location for which the smoothed spline function
#'           is evaluated (0 - 1)
#' @param smooth the smoothing as \code{spar} in \code{\link{smooth.spline}}
#' 
#' @keywords utilities
#' @export
.getval <- function(ind, val, predind, smooth) {
    # Function to get smoothed values in the interval 0 - 1    
    fitfun  <- smooth.spline(ind, val, spar=smooth)
    fitval  <- predict(fitfun, predind)$y
    if (min(fitval) < 0) fitval <- fitval - min(fitval)
    if (max(fitval) > 1) fitval <- fitval / max(fitval)
    fitval
}

#' @rdname colfun
#' @param x Number of colors 
#' @param tmp Array of rgb values (dim=c(ncolors,3)) along which
#'            the color sequence is interpolated
#' @param start,stop the range for which colours should be interpolated (0 - 1)
#' @param log logical, should interval between \code{start} and \code{stop}
#'            be divided linearly?
#' @export
.colseq <- function(x, tmp, start=0, stop=1, log=F, smooth=0.4){
    # .colseq computes a color sequence of x rgb colors along the indicated
    # colors by tmp (Attention: 'ind' is set to .spacing to match the 
    # colours provided by .bluewhite, etc.)
    #
    # Arguments:
    # x         Number of colors 
    # tmp       Array of rgb values (dim=c(ncolors,3)) along which
    #           the color sequence is interpolated
    # start/stop the range for which colours should be interpolated (0 - 1)
    # log       logical, defaults to FALSE, interval between start
    #           and stop is divided linearly
    # smooth    smoothing for smooth.spline in .getval
    #
    # Author    jonas.bhend -at- gkss.de
    #

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

#' @rdname colfun
#' @param n number of colours to return
#' @export
rbfun <- function(n, start=0.4, log=F){
    # rbfun computes a symmetric redblue color scale with x colors
    numx <- floor(n/2)
    if (2*numx == n){
        colour <- c(.colseq((numx+1), .bluewhite, start=start, log=log)[1:numx], 
                    rev(.colseq((numx+1), .redwhite, start=start, log=log)[1:numx]))
    } else {
        colour <- c(.colseq((numx+1), .bluewhite, start=start, log=log), 
                    rev(.colseq((numx+1), .redwhite, start=start, log=log)[1:numx]))
    }
    colour
}

#' @rdname colfun
#' @export
gbfun <- function(n, start=0.4, log=F){
    # rbfun computes a symmetric greenbrown color scale with x colors
    numx <- floor(n/2)
    if (2*numx == n){
        colour <- c(.colseq((numx+1), .brownwhite, start=start, log=log)[1:numx], 
                    rev(.colseq((numx+1), .greenwhite, start=start, log=log)[1:numx]))
    } else {
        colour <- c(.colseq((numx+1), .brownwhite, start=start, log=log), 
                    rev(.colseq((numx+1), .greenwhite, start=start, log=log)[1:numx]))
    }
    colour
}

#' @rdname colfun
#' @export
bbfun <- function(n, start=0.4, log=F){
    # rbfun computes a symmetric bluebrown color scale with x colors
    numx <- floor(n/2)
    if (2*numx == n){
        colour <- c(.colseq((numx+1), .brownwhite, start=start, log=log)[1:numx], 
                    rev(.colseq((numx+1), .bluewhite, start=start, log=log)[1:numx]))
    } else {
        colour <- c(.colseq((numx+1), .brownwhite, start=start, log=log), 
                    rev(.colseq((numx+1), .bluewhite, start=start, log=log)[1:numx]))
    }
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
#' @aliases .soil .water .hsurf .gpcc .redwhite .bluewhite .brownwhite
#'          .greenwhite
#' @name geocolours
NULL
