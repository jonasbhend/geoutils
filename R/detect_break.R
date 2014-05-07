#' Detect step changes
#' 
#' This function checks whether there are step changes present in a time series
#' at a pre-specified time (generally 2005 for CMIP5 historical to RCP scenario
#' transition).
#' 
#' @param x input object of type NetCDF
#' @param breaktime last time step (year) of first period
#' @param n number of years in climatology (before and after break point)
#' @param detrend logical, should linear trend be accounted for?
#' 
#' @keywords utilities
#' @export
detect_break <- function(x, breaktime=2005, n=50, detrend=TRUE){
  xtime <- attr(x, 'time')
  ## rudimentary checking
  if (is.null(xtime))    stop('Break detection on this object not possible (no time attribute)')
  if (breaktime > max(xtime - n) | breaktime < min(xtime + n - 1)){
    stop('Insufficient time steps to perform break detection')
  }
  ## select the two periods
  ndims <- length(dim(x))
  ## compute the standardised difference
  if (detrend){
    dd <- select_time(x, breaktime - n + 1, breaktime+n)
    z <- apply.NetCDF(dd, 1:(ndims-1), function(x) {
      if (any(!is.na(x))){
        xmn1 <- predict(lm(x[1:n] ~ seq(1,n)), interval='conf', level=pnorm(1))
        xmn2 <- predict(lm(x[n + 1:n] ~ seq(1,n)), interval='conf', level=pnorm(1))
        xsd <- sqrt((diff(xmn1[nrow(xmn1), 2:1])**2 + diff(xmn2[1,2:1])**2) / 2)
        xsd[xsd == 0] <- 1
        out <- (xmn1[nrow(xmn1),1] - xmn2[1,1]) / xsd
        return(out)        
      } else {
        return(NA)
      }
    })
  } else {
    ## compute standardised differences
    d1 <- select_time(x, breaktime - n + 1, breaktime)
    d2 <- select_time(x, breaktime + 1, breaktime + n)
    z <- (time_average(d2, n, type='start') - time_average(d1, n, type='start')) / as.vector(sqrt((apply.NetCDF(d1, 1:(ndims - 1), sd, na.rm=T)**2 + apply.NetCDF(d2, 1:(ndims - 1), sd, na.rm=T)**2) / n / 2))
  }
  ## compute probability of these standardised differences occuring (assuming iid)
  zprob <- dt(z, df=2*n - 2, log=T)
  ## combine p-values using Fisher's method
  nind <- 1:(ndims - 2) ## assuming last is seasons and second last is regions
  zfish <- -2 * apply(zprob, nind, sum, na.rm=T)
  ## number of non-NULL t-tests
  kfish <- apply(!is.na(zprob), nind, sum, na.rm=T)
  zfchi <- dchisq(zfish, df=kfish, log=TRUE)
  out <- list(pval=zprob, fisher=zfchi)
  return(out)
}