#' Merge objects of type NetCDF along the time dimension
#' 
#' This function provides support to merge objects of type NetCDF.
#' If overlapping time steps are present, the time step from the
#' second argument are used. If there is a gap between the objects,
#' missing values will be infilled with regularly spaced time steps 
#' in the gap.
#' 
#' @param x input object of class 'NetCDF' 
#' @param y input object of class 'NetCDF'
#' 
#' @keywords utils
#' @export

merge.NetCDF <- function(x,y=NULL){
  xtime <- attr(x, 'time')
  ytime <- attr(y, 'time')
  nseas <- median(c(table(xtime %/% 1), table(ytime %/% 1)))
  atime <- seq(min(c(xtime, ytime)), max(c(xtime, ytime)), by=1/nseas)
  ## relax on constraint for different number of runs
  if (!is.null(rownames(x))){
    unirownames <- rownames(x)[rownames(x) %in% rownames(y)]
    if (length(unirownames) > 0){
      xx <- collapse2mat(x, first=TRUE)
      rownames(xx) <- rownames(x)
      yy <- collapse2mat(y, first=TRUE)
      rownames(yy) <- rownames(y)
      odims <- dim(x)
      odims[1] <- length(unirownames)
      odims[length(odims)] <- length(atime)
      out <- collapse2mat(array(NA, odims), first=FALSE)
      out[,atime %in% xtime] <- as.vector(xx[unirownames,])
      out[,atime %in% ytime] <- as.vector(yy[unirownames,])
      out <- array(out, odims)
      rownames(out) <- unirownames
    } else {
      out <- NULL
    }
  } else if (all(dim(x)[-length(dim(x))] == dim(y)[-length(dim(y))])){
    if (nrow(x) > 0){
      out <- array(NA, c(dim(x)[-length(dim(x))], length(atime)))
      odims <- dim(out)
      out <- collapse2mat(out, first=FALSE)
      out[,atime %in% xtime] <- as.vector(x)
      out[,atime %in% ytime] <- as.vector(y)
      out <- array(out, odims)
    } else {
      out <- NULL
    }
  } else {
    stop('Dimension missmatch when trying to merge netcdf output')
  }
  if (!is.null(out)){
    attnames <- setdiff(names(attributes(x)), c('dim', 'dimnames', 'time'))
    for (atn in attnames){
      attr(out, atn) <- attr(x, atn)
    }
    attr(out, 'time') <- atime
  } 
  return(out)
}
