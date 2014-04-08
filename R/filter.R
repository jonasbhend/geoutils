#' @name fullfilter
#' @aliases applyfilter
#' @aliases get_response
#' @aliases desired_response
#' @aliases get_filter_weights
#' @aliases make_filter
#' 
#' @title Filter functions for 'NetCDF' objects
#'
#' @description
#' \code{fullfilter} applies a moving average filter to the full time
#' series (no missing values, decreasing filter length at start and end).
#' \code{applyfilter} applies \code{fullfilter} to the time dimension
#' of an object of class 'NetCDF'. \code{make_filter} helps you set up a
#' low- band- or highpass filter with the desired response specified in
#' \code{desired_response}. \code{get_filter_weights} provides the weights
#' given a cutoff frequency. \code{get_response} can be used to illustrate
#' and analyse the filter response.
#'
#' @param x time series or 'NetCDF' object
#' @param h filter weights
#' @keywords utilities
#' 
#' @rdname fullfilter
#' @export
fullfilter <- function(x, h){
  if (is.vector(x)){
    
    xout <- filter(x, h)
    nh <- length(h)
    istart <- -(nh%/% 2)
    iend <- nh%/%2 - (nh - 1)%%2
    sumh <- sum(h)
    for (i in which(is.na(xout))){
      ind <- i + istart:iend
      ii <- which(ind > 0 & ind <= length(xout))
      ii <- ii[which(!is.na(x[ind[ii]]))]
      xout[i] <- sum(x[ind[ii]]*h[ii])/sum(h[ii])*sumh
    }
    
  } else {
    stop('filter for non-vectors not implemented')
  }
  xout
}

#' @rdname fullfilter
#' @export
applyfilter <- function(x, h){
  dims <- dim(x)
  ntime <- dims[length(dims)]
  nseas <- median(table(floor(attr(x, 'time'))))
  ##if (length(dim(x)) != 2) stop('Not implemented yet')
  xtmp <- collapse2mat(array(x, c(length(x)/ntime, nseas, ntime/nseas)), first=FALSE)
  xout <- array(NA, dim(xtmp))
  ## run analysis only for non-empty time series
  xout[apply(!is.na(xtmp), 1, any),] <- t(apply(xtmp[apply(!is.na(xtmp), 1, any),,drop=F], 1, fullfilter, h))
  ## change attributes (including dimensions)
  attributes(xout) <- attributes(x)
  return(xout)
}

#' @rdname fullfilter
#' @param accuracy numerical accuracy to which the response should be computed
#' @export
get_response <- function(h, accuracy=1e-5){
  # h contains filter weights
  
  nh  <- length(h)
  if (nh %% 2 == 0) stop("uneven number required")
  a0  <- h[ceiling(nh/2)]
  a   <- h[ceiling(nh/2) + 1:floor(nh/2)]
  k   <- 1:length(a)
  
  # frequencies
  w   <- seq(0,0.5,accuracy)
  cos.arg <- t(t(2*pi*k)) %*% w    
  cw  <- a0 + 2* apply(a * cos(cos.arg), 2, sum)
  
  out <- cbind(w,cw)
  colnames(out) <- c("1/time step", "response")
  out
}


#' @rdname fullfilter
#' @param fr cutoff frequencies
#' @param type one of lowpass, bandpass, or highpass  
#' @export
desired_response <- function(fr, type="lowpass", accuracy=1e-5){
  w   <- seq(-0.5,0.5,accuracy)
  res <- rep(0, length(w))
  
  if (type == "lowpass"){
    res[abs(w) < 1/fr] <- 1
  } else if (type == "highpass"){
    res[abs(w) > 1/fr] <- 1
  } else if (type == "bandpass"){
    fr  <- sort(1/fr)
    res[abs(w) > fr[1] & abs(w) < fr[2]] <- 1
  }
  
  out <- cbind(w, res)
  colnames(out) <- c("1/time step", "response")
  out
}

#' @rdname fullfilter
#' @param n number of time steps used in digital filter
#' @export
get_filter_weights <- function(n, fr, accuracy=1e-5){
  
  tmp <- desired_response(fr, type='lowpass', accuracy=accuracy)
  w   <- tmp[,1]
  cw  <- tmp[,2]        
    
  ak  <- rep(0, n+1)
  hk  <- rep(0, n+1)
  for (i in 0:n){
    ak[i+1] <- sum(cw*cos(2*pi*i*w)*accuracy)
    hk[i+1] <- 0.5*(1+cos(pi*i/(n+1)))
  }
  
  h   <- c(rev(hk), hk[1:n+1])
  a   <- c(rev(ak), ak[1:n+1])
  
  ha  <- h*a
  #ha  <- a
  ha  <- ha/sum(ha)
  ha
}

#' @rdname fullfilter
#' @param ... additional arguments passed to \code{get_filter_weights}
#' @examples
#' xx <- rnorm(100) + sin(seq(0,4*pi, length=100))
#' plot(xx, type='l')
#' hh <- make_filter(31, 10, type='lowpass')
#' lines(filter(xx, hh), lwd=2)
#' plot(get_response(hh), type='l')
#' abline(v=1/10, lty=2, lwd=2)
#' @export
make_filter <- function(n, fr, type="lowpass", ...) {
  
  
  fr  <- sort(fr)
  
  if (n %% 2 == 0) stop("Uneven number of coefficients required")
  u   <- rep(0,n)
  u[ceiling(n/2)] <- 1
  a   <- get_filter_weights( n=(n-1)/2, fr=fr[1], ...)
  
  if (type == "highpass"){
    a   <- u - a
  } else if (type == "bandpass"){
    a   <- a - get_filter_weights( n=(n-1)/2, fr=fr[2], ...)
  }
  
  as.vector(a)
}


