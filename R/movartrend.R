#' @name movartrend
#' @aliases trendresid trendcov artrend convertx
#' 
#' @title
#' Compute trends on NetCDF objects
#' 
#' @description
#' Compute trends on NetCDF objects where \code{movartrend} computes moving 
#' linear regressions on the time series, \code{trendresid} and \code{trendcov} 
#' compute simple linear regression with the latter returning the covariance
#' matrix of the standard errors (i.e. for multivariate analyses).
#' 
#' The function \code{convertx} is used to expand an object of type NetCDF
#' to fit another object of type NetCDF. The standard use case would be to 
#' expand a NetCDF object with ensemble mean global mean temperature to an
#' object with multiple models and simulations. \code{convertx} is used in 
#' \code{movartrend} to convert a potential predictor for the regression.
#' 
#' 
#' @param data object of type NetCDF containing the time series
#' @param trndln length of trend to be computed (in years if time 
#' attribute is present)
#' @param x optional covariate to be used as predictor
#' @param missthresh threshold to indicate what fraction of non-missing values is 
#' at least required to run the regression
#' @param robust logical, should robust regression be used?
#' @param gaps logical, should gaps be allowed for?
#' @param plotable logical, should output be rearrange to make it plotable?
#' @param resid logical, should residuals of regression be returned?
#' 
#' @examples
#' ## result should be roughly 1
#' print(movartrend(t(1:100 + rnorm(100, sd=10)), trndln=100)$trend[1])
#' 
#' @useDynLib geoutils
#' 
#' @keywords utilities
#' @export
movartrend <- function(data, trndln, x=NULL, missthresh = 0.8, robust = FALSE, gaps = FALSE, plotable=TRUE, resid=FALSE){
  
  if (length(dim(data)) > 1){
    dims <- dim(data)
    dat  <- collapse2mat(data, first=FALSE)
    nseas <- 1
    if ('time' %in% names(attributes(data))){
      tim <- attr(data, 'time')
      nseas <- median(table(floor(tim)))
      dat <- collapse2mat(array(dat, c(nrow(dat), nseas, floor(ncol(dat)/nseas))), first=FALSE)
    }
    missval <- as.double(-1e20)
    dat[is.na(dat)] <- missval
    
    y <- array(as.double(dat), dim(dat))
    y.out <- array(missval, dim(dat))
    se0 <- array(missval, dim(dat))
    se1 <- array(missval, dim(dat))
    se2 <- array(missval, dim(dat))
    trndln <- as.integer(trndln)
    missval <- as.double(-1e20)
    n <- as.integer(dim(dat))
    missthresh <- as.double(missthresh)
    robust <- as.integer(robust)
    gaps <- as.integer(gaps)
    x <- convertx(x, data=data, dat=dat)
    dat[is.na(x)] <- missval
    x[is.na(x)] <- missval
    
    tmp <- .Fortran("artrend",
                    n = n,
                    y = y,
                    x = x, 
                    trndln = trndln, 
                    yout = y.out,
                    yresid = y.out,
                    seout = y.out,
                    se0 = se0,
                    missval = missval, 
                    missthresh = missthresh, 
                    robust = robust, 
                    gaps = gaps)
    tmp <- lapply(tmp, function(x) {x[x == missval] <- NA ; return(x)})
    tnames <- names(tmp)
    tind <- which(tnames %in% c('yout', 'seout', 'se0'))
    if (resid & ncol(tmp$yout) == trndln) tind <- c(tind, which(tnames == 'yresid'))
    tmp <- tmp[tind]
    if (resid & ncol(tmp$yout) == trndln) {
      names(tmp) <- c('trend', 'se', 'se0', 'resid')
    } else {
      names(tmp) <- c('trend', 'se', 'se0')
    }
    if (ncol(tmp$trend) == trndln){
      tmpout <- lapply(tmp, function(x) x[,trndln])
      if (any(names(tmpout) == 'resid')) tmpout$resid <- tmp$resid
      tmp <- tmpout
      dims[length(dims)] <- nseas
    } else if (!plotable) {
      dims <- c(dims[-length(dims)], nseas, ncol(tmp$trend))
    }
    
    out <- lapply(tmp, function(x) {
      if (prod(dims) == length(x)){
        x <- array(x, dims)
      } else {
        x <- array(x, c(dims[-length(dims)], nseas*trndln))
      }
      for (atn in setdiff(names(attributes(data)), 'dim')){
        attr(x, atn) <- attr(data, atn)
      }
      if ('time' %in% names(attributes(data))){
        if (length(attr(data, 'time')) == dim(x)[length(dim(x))]){
          attr(x, 'time') <- attr(data, 'time')
        } else {
          attr(x, 'time') <- attr(data, 'time')[1:dim(x)[length(dim(x))]]
        }
      }
      return(x)
    })
    
    return(out)
  }
  
}

#'@rdname movartrend
#'@export
trendcov <- function(data, x=NULL, missthresh = 0.8){

  if (length(dim(data)) > 1){
    dims <- dim(data)
    dat  <- collapse2mat(data, first=FALSE)
    nseas <- 1
    if ('time' %in% names(attributes(data))){
      tim <- attr(data, 'time')
      nseas <- median(table(floor(tim)))
      dat <- collapse2mat(array(dat, c(nrow(dat), nseas, ncol(dat)/nseas)), first=FALSE)
    }
    dims <- c(dims[-length(dims)], nseas)
    missval <- as.double(-1e20)
    x <- convertx(x, data=data, dat=dat)
    dat[is.na(dat)] <- missval
    x[is.na(x)] <- missval
    tmp <- .Fortran("trendcov",
                    n = as.integer(dim(dat)),
                    y = array(as.double(dat), dim(dat)),
                    x = array(as.double(x), dim(x)), 
                    beta = as.double(rep(0, nrow(dat))),
                    betacov = array(as.double(0), rep(nrow(dat), 2)),
                    missval = as.double(missval), 
                    missthresh = as.double(missthresh)
    )
    
    out <- lapply(tmp[names(tmp) %in% c('beta', 'betacov')], function(x) {x[x == missval] <- NA ; return(x)})
    
    out$beta <- array(out$beta, dims) 
    for (atn in setdiff(names(attributes(data)), 'dim')){
      attr(out$beta, atn) <- attr(data, atn)
    }
    if ('time' %in% names(attributes(data))){
      attr(out$beta, 'time') <- range(floor(attr(data, 'time')))
    }
    return(out)
  }
  
}

#'@rdname movartrend
#'@export
trendresid <- function(data, x=NULL, missthresh = 0.8){

  if (length(dim(data)) > 1){
    dims <- dim(data)
    dat  <- collapse2mat(data, first=FALSE)
    nseas <- 1
    if ('time' %in% names(attributes(data))){
      tim <- attr(data, 'time')
      nseas <- median(table(floor(tim)))
      dat <- collapse2mat(array(dat, c(nrow(dat), nseas, ncol(dat)/nseas)), first=FALSE)
    }
    dims <- c(dims[-length(dims)], nseas)
    missval <- as.double(-1e20)
    x <- convertx(x, data=data, dat=dat)
    dat[is.na(dat)] <- missval
    x[is.na(x)] <- missval
    tmp <- .Fortran("trendresid",
                    n = as.integer(dim(dat)),
                    y = array(as.double(dat), dim(dat)),
                    x = array(as.double(x), dim(x)), 
                    beta = as.double(rep(0, nrow(dat))),
                    dw = as.double(rep(0, nrow(dat))),
                    se0 = as.double(rep(0, nrow(dat))),
                    seout = as.double(rep(0, nrow(dat))),
                    missval = as.double(missval), 
                    missthresh = as.double(missthresh))
    
    out <- lapply(tmp[names(tmp) %in% c('beta', 'y', 'dw', 'se0', 'seout')], function(x) {x[x == missval] <- NA ; return(x)})
    
    out$yresid <- array(out$y, dim(data))
    attributes(out$yresid) <- attributes(data)
    out$y <- NULL
    for (n in setdiff(names(out), 'yresid')){
      out[[n]] <- array(out[[n]], dims) 
      for (atn in setdiff(names(attributes(data)), 'dim')){
        attr(out[[n]], atn) <- attr(data, atn)
      }
      if ('time' %in% names(attributes(data))){
        attr(out[[n]], 'time') <- range(floor(attr(data, 'time')))
      }
    }
    return(out)
  }
}    

#' @rdname movartrend
#' @export
artrend <- function(data, x=NULL, missthresh=0.8){
  if (sum(!is.na(data)) > 5 & mean(!is.na(data)) > missthresh){
    N <- length(data)
    if (is.null(x)) x <- seq(1,N)
    if (any(is.na(data))){
      mask <- !is.na(data)
      x <- x[mask]
      data <- data[mask]
      N <- sum(mask)
      warning('Missing values removed, computation of AR coefs not reliable')
    }
    tlm <- lm(data ~ x)
    out <- summary(tlm)$coef[2,]
    ac  <- acf(tlm$res, plot=F)$acf
    neff1 <- nefffun(N, acf=ac[2])
    neff2 <- nefffun(N, acf=ac[2:3])
    out <- summary(tlm)$coef[2,]
    out[3] <- sqrt(out[2]**2 * (N-tlm$rank)/(neff1 - tlm$rank))
    out[4] <- sqrt(out[2]**2 * (N-tlm$rank)/(neff2 - tlm$rank))
  } else {
    out <- rep(NA, 4)
  }
  names(out) <- c('Estimate', 'SE', 'SE - AR(1)', 'SE - AR(2)')
  return(out)
}

#' @rdname movartrend
#' @param x object of type NetCDF to be expanded to fit
#' dimensions of \code{dat}
#' @param data object of type NetCDF to be used to match
#' attributes (e.g. time steps etc.)
#' @param dat object to infer dimensions of output
#' @export
convertx <- function(x=NULL, data, dat=data){
  if (is.null(x)){
    if ('time' %in% names(attributes(data))){
      tim <- attr(data, 'time')
      nspace <- length(data)/length(tim)
      x <- array(rep(tim, each=nspace), dim(dat))
    } else {
      x <- array(rep(1:ncol(dat), each=nrow(dat)), dim(dat))
    }
  } else {
    if (length(data) %% length(x) != 0) stop("Dimensions of x and data don't match")
    if (length(x) == length(data)){
      x <- array(x, dim(dat))
    } else if (length(x) == dim(data)[length(dim(data))]){
      x <- array(rep(x, each=prod(dim(data)[-length(dim(data))])), dim(dat))
    } else {
      ## find the correct dimensions to expand in data
      xdims <- dim(x)
      ddims <- dim(data)
      if (length(xdims) == length(ddims) - 1){
        xi <- which(ddims %in% xdims)
        if (length(xi) == length(ddims)) {
          warning('Ambiguity in dimensions --> assuming first has to be expanded')
          xi <- xi[2:length(xi)]
        }
        xout <- aperm(array(rep(x, length=length(data)), c(ddims[xi], ddims[-xi])), c(xi, setdiff(seq(along=ddims), xi)))
        x <- array(xout, dim(dat))
      } else {
        stop(paste('Case not supported: dimensions of data', ddims, 'dimension of x', xdims))
      }
    }
  }
  ## rewrite the attributes
  attributes(x) <- attributes(dat)
  return(x)
}

