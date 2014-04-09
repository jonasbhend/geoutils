#' Plot shaded area between lines
#' 
#' Function to simplify plotting shaded areas between two lines.
#' Functionality to deal with missing values is added (breaks in
#' the shaded areas where there are missing values).
#' 
#' @param y1,y2 vectors for vertical lines
#' @param x vector with x intersects (optional)
#' @param horizontal logical, x-y axis inverted if \code{FALSE}
#' @param ... additional arguments passed to \code{\link{polygon}}
#' 
#' @keywords plot
#' @export
vert.poly <- function(y1, y2, x=seq(along=y1), horizontal=TRUE, ...) {
  # plots polygons using polygon having lower bounds y1
  # and upper bounds y2 along x
  
  y1[abs(y1) == Inf]  <- NA
  y2[abs(y2) == Inf]  <- NA
  
  ndim <- length(y1)
  if (ndim == length(y2) & ndim == length(x)){
    
    mask    <- which(!is.na(y1) & !is.na(y2))
    nmask   <- length(mask)
    
    # find subsets for which to draw a polygon
    gaps    <- c(0,which(diff(mask) > 1),nmask)
    ngaps   <- length(gaps)-1
    
    for (i in 1:ngaps){
      index   <- mask[gaps[i]+1]:mask[gaps[i+1]]
      y.tmp   <- c(y1[index],y2[rev(index)])
      x.tmp   <- c(x[index], x[rev(index)])
      if (horizontal){
        polygon(x.tmp, y.tmp, ...)
      } else {
        polygon(y.tmp, x.tmp, ...)
      }
    }
    
  } else {
    print("dimensions of y1, y2, and x do not match")
  }
}    
