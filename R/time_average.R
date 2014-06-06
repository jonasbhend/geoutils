#' Create n-year averages of an NetCDF object
#'
#' This function computes n-year averages from an existing NetCDF time series 
#' read in with \code{\link{readNetCDF}} 
#'
#' @param x array of class 'NetCDF'
#' @param n number of time steps to average over
#' @param offset offset from which to start (in years)
#' @param type Alternative to specifying offset, 
#' one of start (the default), even or end, see details
#' @details
#' The automated way to set the offset works as following.
#' \code{start} starts the \code{n}-year averages at the beginning
#' of the time series, \code{end} sets the offset so that the last
#' average ends in the last year of the time series, \code{even} 
#' sets the offset so that \code{modulo(year, n) == 0} for the first year
#' in the first average (i.e. for \code{n=5} averages such as 1910-1914).
#' 
#' The time attribute reflects the average over the same block as the average
#' of the data, that is \code{n}-yearly averages will be centered.
#' 
#' @examples
#' tas <- readNetCDF(system.file("extdata", "annual_CRUTEMv3_1961-90.nc", package="geoutils"), varname="temp")
#' 
#' ## compute 5-yearly averages
#' tas.avg <- time_average(tas, 5, type='start')
#' 
#' ## get a grid point with little missing values
#' si <- which.max(apply(!is.na(tas), 1, sum))
#' 
#' ## plot the data
#' plot(tas, type='ts', si=si, lwd=1)
#' plot(tas.avg, type='ts', si=si, col=2, lwd=2, add=T)
#' 
#' @keywords utilities
#' @export

time_average <- function(x, n, offset=NULL, type='start'){
  
  if (!is.null(attr(x, "time"))){
    
    # get the time attribute
    tim     <- attr(x, "time")
    nseas   <- sum(tim < tim[1] + 1)
    
    # dimensions in input
    dims    <- dim(x)
    if (is.null(dims)) dims <- length(x)
    
    if (length(tim) != dims[length(dims)]) {
      print(paste("WARNING: time attribute and dimensions mismatch in", get(x)))
    }
    
    # compute offset (offset = 0 for control runs)
    # offset in years
    if (is.null(offset)){
      if (type == 'even') {
        offset  <- (which.max((tim %% n) == 0) - 1)/nseas
      } else if (type == 'start') {
        offset <- 0
      } else if (type == 'end') {
        offset <- (length(tim)/nseas) %% n
      }
    }
    
    # number of timesteps and blocks in output
    ntimes  <- length(tim) - offset*nseas
    ntimes <- floor(ntimes/nseas)
    nblocks <- floor(ntimes/n)
    
    dim.new <- dims
    dim.new[length(dims)] <- nblocks*nseas
    if (dims[length(dims)] == 1){
      warning('Only one time step in input - assuming input is average')
      out <- x
    } else if (length(dim.new) == 1){
      dim.int <- c(1,nseas,floor(dims[length(dims)]/nseas))
      dim.int2<- c(1,nseas, n, nblocks)
      tmp     <- array(x, dim.int)
      out <- array(apply(array(tmp[,,offset+1:(n*nblocks)], dim.int2), c(1,2,4), mean, na.rm=T), dim.new)               
    } else {
      dim.int <- c(prod(dims[1:(length(dims)-1)]), nseas, floor(dims[length(dims)]/nseas))
      dim.int2 <- c(prod(dims[1:(length(dims)-1)]), nseas, n, nblocks)
      tmp     <- array(x, dim.int)
      tmp     <- tmp[,,offset+1:(n*nblocks), drop=F]
      tmp.na  <- (!is.na(tmp))* 1
      tmp[is.na(tmp)] <- 0
      ntmp    <- dim(tmp)[3]
      tmp.out <- tmp[,,seq(1,ntmp, n)]
      na.sum  <- tmp.na[,,seq(1,ntmp,n)]
      for (j in 2:n){
        tmp.out <- tmp.out + tmp[,,seq(j,ntmp,n)]
        na.sum  <- na.sum + tmp.na[,,seq(j,ntmp,n)]
      }
      tmp.out[na.sum != 0]    <- tmp.out[na.sum != 0] / na.sum[na.sum != 0]
      tmp.out[na.sum == 0]    <- NA
      
      out <- array(tmp.out, dim.new)       
    }
    
    attnames <- setdiff(names(attributes(x)), 'dim')
    for (atn in attnames){
      attr(out, atn) <- attr(x, atn)
    }
    # convert the time accordingly
    if (dims[length(dims)] == 1){
      attr(out, 'time') <- attr(x, 'time')
    } else {
      timtmp <- tim[offset*nseas + 1:(n*nblocks*nseas)]
      tim.new <- ceiling(as.vector(apply(array(floor(timtmp), c(nseas, n, nblocks)), c(1,3), mean))) + timtmp[1:nseas] %% 1
      attr(out, "time") <- tim.new      
    }
  } else {
    out <- x
  }
  
  invisible(out)
  
}
