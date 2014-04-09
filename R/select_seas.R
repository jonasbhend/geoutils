#' Select seasons from 'NetCDF' object
#' 
#' Function to select specified season (by fractional year, i.e. 0.5 for
#' summer) from an object of class 'NetCDF'.
#' 
#' @param x object of class 'NetCDF'
#' @param seas fractional year indicating what season to select
#' 
#' @keywords utilities
#' @export
select_seas <- function(x, seas){
  if ('time' %in% names(attributes(x))){
    tim <- attr(x, 'time')
    ntim <- length(tim)
    xdim <- dim(x)
    if (length(xdim) == 0) xdim <- length(x)
    ti <- which(xdim == ntim)
    if (length(ti) > 1) {
      warning('Time dimension probably mixed up')
      ti <- max(ti)
    }
    
    ## find matching times
    tseas <- round(tim %% 1, 5)
    seas <- round(seas, 5)
    timei <- which(tseas == seas)
    if (length(timei) == 0){
      timei <- which((tseas - seas) ** 2 == min((tseas - seas)**2))
      warning('Guessing the right time index')
    }
    if (length(xdim) == 1){
      xout <- x[timei]
    } else {
      perm <- c(ti, setdiff(seq(along=xdim), ti))
      xout <- collapse2mat(aperm(x, perm), first=TRUE)[timei,]
      if (ti == 1){
        reperm <- seq(along=xdim)
      } else {
        reperm <- seq(2, length(xdim)+1)
        reperm[ti] <- 1
        reperm[reperm > ti] <- reperm[reperm > ti] - 1
      }
      outdim <- xdim
      outdim[ti] <- length(timei)
      xout <- aperm(array(xout, outdim[perm]), reperm)
    }
    attnames <- setdiff(names(attributes(x)), 'dim')
    for (atn in attnames){
      attr(xout, atn) <- attr(x, atn)
    }
    if (length(which(tim %% 1 == seas)) > 0){
      attr(xout, 'time') <- tim[timei]
    } else {
      attr(xout, 'time') <- tim[timei] - tim[timei] %% 1 + seas
    }
  } else {
    xout <- x
  }
  return(xout)
}
