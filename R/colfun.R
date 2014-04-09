#' @name colfun
#' @aliases .getval
#' @aliases .colseq
#' @aliases rbfun
#' @aliases gbfun
#' @aliases bbfun
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

.soil   <- c("#E6FAFF","white", "grey", gbfun(6)[c(3:1,4:6)], .water)

.data <- t(array(c(
  0, 0, 0,
  2, 31, 64,
  5, 48, 97,
  23, 62, 132,
  33, 102, 172,
  67, 147, 195,
  146, 197, 222,
  209, 229, 240,
  253, 253, 253,
  253, 219, 199,
  244, 165, 130,
  214, 96, 77,
  178, 24, 43,
  138, 14, 33,
  103, 0, 31,
  68, 0, 28,
  0, 0, 0),
  dim=c(3,17)))/255

.hsurf  <- t(array(c(
  60,  179,  113,
  157,  210,  156,
  207,  229,  174,
  242,  244,  193,
  230,  226,  148,
  216,  211,  114,
  165,  152,   72,
  121,   79,   31,
  81,   63,   22),
  dim=c(3,9)))/255

.gpcc   <- t(array(c(
  255, 252, 127,
  255, 250, 0,
  132, 250, 127,
  22, 247, 0,
  18, 194, 0,
  15, 165, 0,
  14, 162, 255,
  9, 109, 255,
  0, 39, 255,
  139, 50, 222,
  159, 31, 161,
  189, 37, 191,
  218, 43, 221,
  253, 0, 255),
  dim=c(3,14)))/255

# define the standard colour bars
.bluewhite  <- .data[1:9,]
.redwhite   <- .data[17:9,]
.greenwhite <- .data[17:9,c(2,1,3)]
.brownwhite <- .data[1:9,c(3,2,1)]
.water        <- "#ABE1FA"

rm(.data)
