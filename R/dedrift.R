#' Function to remove drift from model runs
#' 
#' Simple drift correction of model runs given a corresponding control run
#' 
#' @param data (list of) input time series of type NetCDF
#' @param ctl (list of) control run time series of type NetCDF
#' 
#' @examples
#' ## synthetic historical and control time series
#' xx <- t(rnorm(400) + c(-4,0,3,1) + seq(0,sqrt(10), length=400)**2)
#' attr(xx, 'time') <- seq(1900, by=0.25, length=400)
#' class(xx) <- 'NetCDF'
#' 
#' ## synthetic control runs
#' xctl <- t(as.vector(outer(c(-5,-3,-2,-1), log(seq(1, 100, length=300)), '*')) + rnorm(1200))
#' attr(xctl, 'time') <- seq(1750, by=0.25, length=length(xctl))
#' class(xctl) <- 'NetCDF'
#' 
#' ## without dedrifting
#' plot(xx, type='ts', seas=T, main='Synthetic time series (forced simulation)')
#' 
#' ## control run
#' plot(xctl, type='ts', seas=T, main='Synthetic control')
#' 
#' ## dedrifted forced run
#' plot(dedrift(xx, xctl), type='ts', seas=T, main='Dedrifted forced simulation') 
#' plot(xx, type='ts', seas=T, lwd=1, lty=2, add=T)
#' 
#' @keywords utilities
#' @export
dedrift <- function(data, ctl){
  if (is.list(data)){
    vars <- names(data)
    out <- list()
  } else {
    vars <- 1
  }
  for (var in vars){
    if (is.list(data)){
      dat <- data[[var]]
      ct <- ctl[[var]]
    } else {
      dat <- data
      ct <- ctl
    }
    ## select overlapping time periods
    ctime <- attr(ct, 'time')
    dtime <- attr(dat, 'time')
    ctim <- range(floor(ctime))
    dtim <- range(floor(dtime))
    if (all(dtime %in% ctime)){
      ct <- select_time(ct, dtim[1], dtim[2])
    } else {
      ct <- select_time(ct, ctim[1], ctim[1] + diff(dtim))
    }
    ## make sure the data have the same units
    ct <- convertx(ct, dat)
    ## compute the trends from the control
    ct.trend <- movartrend(ct, diff(dtim)+1)$trend
    ## expand to get residuals for individual time steps
    xtim <- dtim[1]:dtim[2] - mean(dtim)
    ct.resid <- as.vector(as.vector(ct.trend) * rep(xtim, each=length(ct.trend)))
    ## subtract ct trend
    ddat <- dat - ct.resid
    attr(ddat, 'trend') <- ct.trend
    
    if (is.list(data)){
      out[[var]] <- ddat
    } else {
      out <- ddat
    }
  }
  invisible(out)
}
