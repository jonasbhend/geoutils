#' Compute anomalies with respect to reference period
#' 
#' This function computes anomalies with respect to a given reference period. For
#' the function to work, the reference period must be fully in the time series
#' provided.
#' 
#' @param x object of type NetCDF
#' @param anomaly start and end year of reference period
#' @param relative logical, should relative anomalies (in percent) be computed?
#' @param na.thresh minimum fraction of non-missing values in reference period
#' @param rel.thresh maximum relative anomalies that get computed
#' @param zero.thresh maximum number of zeros in reference for relative anomalies
#' to be computed
#' 
#' @keywords utilities
#' @export
anomalies <- function(x, anomaly=c(1961,1990), relative=FALSE,
                      na.thresh=0.8, rel.thresh=1e6, zero.thresh=1){
  # computes anomalies for the respective time selected
  xtim <- attr(x, 'time')
  if (!is.null(xtim)){
    nseas <- median(table(floor(xtim)))
    xx <- array(x, c(length(x)/length(xtim)*nseas, length(xtim)/nseas))
    if (min(xtim) < anomaly[1] + max(diff(xtim)) & max(xtim) > anomaly[2] - max(diff(xtim))){
      xtmp <- select_time(x, anomaly[1], anomaly[2])
      xttim <- attr(xtmp, 'time')
      if (xttim[1] %% 1 != xtim[1] %% 1){
        xtmp <- array(xtmp, c(length(xtmp)/length(xttim), nseas, length(xttim)/nseas))
        perm <- NULL
        for (i in 1:nseas){
          perm <- c(perm, which(xttim[1:nseas] %% 1 == xtim[i] %% 1))
        }
        xtmp <- array(xtmp[,perm,], c(length(xtmp)/length(xttim)*nseas, length(xttim)/nseas))
      } else {
        xtmp <- array(xtmp, c(length(xtmp)/length(xttim)*nseas, length(xttim)/nseas))
      }
    } else if (min(xtim) < 100){
      xtmp <- xx
    } else {
      stop('Anomaly interval not in data')
    }
    x.mn <- apply(xtmp, 1, mean, na.rm=T)
    x.na <- apply(!is.na(xtmp), 1, sum, na.rm=T)
    x.mn[x.na < ncol(xtmp)*na.thresh] <- NA
    if (relative){
      xout <- (xx - x.mn)/x.mn*100
      xout[apply(xout, 1, max, na.rm=T) > rel.thresh |
             apply(xtmp == 0, 1, sum, na.rm=T)/x.na > zero.thresh,] <- NA
    } else {
      xout <- xx - x.mn
    }
    out <- array(c(xout, rep(NA, length(x) - length(xx))), dim(x))
    for (atn in setdiff(names(attributes(x)), 'dim')){
      attr(out, atn) <- attr(x, atn)
    }
    return(out)
  } else {
    return(x)
  }
}
