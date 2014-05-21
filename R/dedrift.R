#' Function to remove drift from model runs
#' 
#' Simple drift correction of model runs given a corresponding control run
#' 
#' @param data (list of) input time series of type NetCDF
#' @param ctl (list of) control run time series of type NetCDF
#' @param varname variable name (if data is no list)
#' 
#' @keywords utilities
#' @export
dedrift <- function(data, ctl, varname='pr'){
  warning('De-drifting should be used with relative or log-normalized precip data only')
  if (is.list(data)){
    vars <- names(data)
    out <- list()
  } else {
    vars <- varname
  }
  for (var in vars){
    if (is.list(data)){
      dat <- data[[var]]
      ct <- ctl[[var]]
    } else {
      dat <- data
      ct <- ctl
    }
    if (any(dim(dat)[c(1,3)] != dim(ct)[c(1,3)])) stop('Dimensions do not match')
    ct.mn <- apply(ct, c(1,3,4), mean, na.rm=T)
    ct.arr <- array(ct.mn, c(prod(dim(ct.mn)[1:2]), dim(ct.mn)[3]))
    ct.trend <- rep(NA, nrow(ct.arr))
    mask <- apply(!is.na(ct.arr), 1, any)
    ct.trend[mask] <- apply(ct.arr[mask,], 1, function(x) lm(x ~ seq(along=x))$coef[2])
    ct.trend <- array(ct.trend, c(dim(ct.mn)[1], 1, dim(ct.mn)[2]))[,rep(1, dim(dat)[2]),]
    index <- dim(dat)[4]:1
    runs <- apply(apply(!is.na(dat), 1:2, any), 1, sum)
    for (i in 1:dim(dat)[1]){
      for (j in 1:runs[i]){
        for (k in 1:dim(dat)[3]){
          ind <- index
          ind[is.na(dat[i,j,k,])] <- NA
          ind <- ind - mean(ind, na.rm=T)
          dat[i,j,k,] <- dat[i,j,k,] + ind*ct.trend[i,j,k]
        }
      }
    }
    if (is.list(data)){
      out[[var]] <- dat
      attnames <- setdiff(names(attributes(data[[var]])), 'dim')
      for (atn in attnames){
        attr(out[[var]], atn) <- attr(data[[var]], atn)
      }
      attr(out[[var]], 'trend') <- ct.trend[,1,]
    } else {
      out <- dat
      attnames <- setdiff(names(attributes(data)), 'dim')
      for (atn in attnames){
        attr(out, atn) <- attr(data, atn)
      }
      attr(out, 'trend') <- ct.trend[,1,]
    }
  }
  invisible(out)
}
