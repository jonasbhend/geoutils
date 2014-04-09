#' Compute (multi-model) ensemble means
#' 
#' This function computes the (unweighted) ensemble mean for an initial
#' condition ensemble or a multi-model ensemble on the assumption that the
#' input of class \code{NetCDF} is a 3d or 4d array with the first dimension 
#' reflective of models and the second dimension containing runs (in the 4d case).
#' 
#' @param x input array of class 'NetCDF'
#' @param multimodel logical, should first dimension (models) be averaged?
#' 
#' @keywords utilities
#' @examples
#' xtmp <- array(1:5, c(5,3,1,1))
#' class(xtmp) <- 'NetCDF'
#' ensmean(xtmp)[,1,1]
#' ensmean(ensmean(xtmp), multimodel=TRUE)[,1]
#' @export
#' 
ensmean <- function(x, multimodel=FALSE){
  if (is.list(x)){
    out <- lapply(x, function(y) {
      ny <- nrow(y)
      if (ny == 1){
        attr(y, 'nens') <- ny
        return(y)
      } else {
        out <- y[1,,,drop=FALSE]
        for (i in 2:ny) out <- out + as.vector(y[i,,])
        out <- out / ny
        attributes(out) <- c(attributes(out), attributes(y)[! names(attributes(y)) %in% c('dim', 'dimnames')])
        attr(out, 'nens') <- ny
        return(out)
      }})
  } else {
    if (multimodel){
      out <- apply(x, 2:length(dim(x)), mean)
      if (length(dim(x)) == 2) out <- t(out)
    } else {
      out <- apply(x, setdiff(seq(1,length(dim(x))), 2), mean, na.rm=T)
    }
    nens <- apply(apply(!is.na(x), 1:2, any), 1, sum)
    attnames <- setdiff(names(attributes(x)), c('dim', 'dimnames'))
    for (atn in attnames){
      attr(out, atn) <- attr(x, atn)
    }
    attr(out, 'nens') <- nens
  }
  invisible(out)
}  
