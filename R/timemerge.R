#' Merge a list of NetCDF objects along the time dimension
#' 
#' Function to merge different NetCDF objects (e.g. historical and scenario simulation)
#' 
#' @param x list of NetCDF objects to be merged
#' 
#' @keywords utilities
#' @export
timemerge <- function(x) {
  if (!is.list(x)) stop('timemerge works on list objects only')
  times <- sort(unique(unlist(lapply(x, attr, 'time')))) 
  times <- seq(min(times), max(times), by=1/median(table(times %/% 1)))
  xdims <- sapply(x, function(y) dim(y))
  if (!all(xdims[seq(1, nrow(xdims) - 1),] == xdims[seq(1, nrow(xdims) - 1),1])) stop("Non-time dimensions in timemerge don\'t match")
  out <- array(NA, c(prod(xdims[seq(1, nrow(xdims) - 1),1]), length(times)))
  for (i in seq(along=x)){
    out[,which(times %in% attr(x[[i]], 'time'))] <- collapse2mat(x[[i]], first=FALSE)[,which(attr(x[[i]], 'time') %in% times)]
  }
  out <- array(out, c(xdims[seq(1, nrow(xdims) - 1),1], length(times)))
  attnames <- setdiff(names(attributes(x[[1]])), 'dim')
  for (atn in attnames){
    attr(out, atn) <- attr(x[[1]], atn)
  }
  attr(out, 'time') <- times
  return(out)
}

