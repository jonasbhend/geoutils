#' Plot a colourbar
#' 
#' Function to plot a colourbar with specified colours and leves
#' or with output from \code{plot.NetCDF}. 
#' 
#' @param levs vector of shading boundaries or list containing \code{levs}, 
#'            \code{cols}, \code{units} and \code{sea.col}
#' @param cols vector of colours
#' @param units character string describing units (optional)
#' @param side which side should axis labels be drawn on? \code{side = 2,4} 
#'             will produce a vertically aligned colour bar        
#' @param xlab,ylab axes labels label (alternative to describe the units)
#' @param labels alternative labels to be used instead of levs
#' @param nn spacing of labels to be drawn, \code{nn = 2} draws every other label only
#' @param center logical, should labels refer to colour mid-points?
#' @param cex.axis character expansion for axis labelling
#' @param sea.col colour for missing values, masked values
#' @param sea.name description of missing/masked values
#' 
#' @keywords utilities
#' @export
#' 
plot_colourbar <- function(levs, cols=NULL, units=NULL, side=1, ylab="", labels=NULL, xlab="", nn=1, center=F, cex.axis=par('cex.axis'), sea.col=NULL, sea.name='N/A', ...){
  
  ## if the first argument is a list, extract the elements
  if (is.list(levs)){
    if (is.null(cols)) cols <- levs$col
    if (is.null(units)) units <- levs$unit
    if (is.null(sea.col)) sea.col <- levs$sea.col
    levs <- levs$lev
  }

  ## colour for mask not part of the standard colour ramp
  sea.add     <- 0
  if (!is.null(sea.col)) {
    cols    <- c(sea.col, cols)
    sea.add <- 1
  }
  ncols       <- length(cols)
  lev.arr     <- array(1:ncols, c(1, ncols))
  if (side %in% c(1,3)){
    lev.arr <- t(lev.arr)
  }    
  
  if (side %in% c(1,3)){
    image(1:ncols, 1, lev.arr, axes=F, breaks=1:(ncols+1)-0.5, col=cols,
          ylab=ylab, xlab=xlab, ...)
    abline(v=1:(ncols-1) + 0.5)
  } else {
    image(1, 1:ncols, lev.arr, axes=F, breaks=1:(ncols+1)-0.5, col=cols,
          ylab=ylab, xlab=xlab, ...)
    abline(h=1:(ncols-1) + 0.5)
  }
  if (center){
    at.lev  <- seq(1+sea.add,ncols,nn)
    if (is.null(labels)){
      labels <- ((levs[seq(1,ncols-sea.add)] + levs[seq(2,ncols+1-sea.add)])/2)[at.lev - sea.add]
    }
    axis(side, at=at.lev, labels=labels, las=1, tick=F, cex.axis=cex.axis)
  } else {
    at.lev  <- seq(2+sea.add,ncols,nn)
    if (is.null(labels)){
      labels  <- levs[at.lev - sea.add]
    }
    axis(side, at=at.lev-0.5,labels=labels, las=1, cex.axis=cex.axis) 
  }
  if (!is.null(sea.col)){
    axis(side, at=1, labels=sea.name, las=1, tick=F, cex.axis=cex.axis)
  }
  box()
  
  # plot unit-string in axis
  if (!is.null(units)){
    if (strwidth(units) > 0.9){
      axis(side, at=par('usr')[2]-0.5, labels=units, 
           tick=F, las=1, cex.axis=cex.axis, hadj=0)
    } else {
      axis(side, at=par('usr')[2], labels=units, 
           tick=F, las=1, cex.axis=cex.axis)
    }
  }
}


