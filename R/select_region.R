#' Select region(s)
#' 
#' Function to select spatial entities (2nd last dimension in NetCDF array).
#' 
#' @param x input array of class NetCDF
#' @param regi index or name of region(s) to select
#' 
#' @details
#' Named regions to select only work if the NetCDF object contains an attribute
#' 'regions' that contains the region names referred to in \code{regi}.
#' 
#' @keywords utilities
#' @export
select_region <- function(x, regi=1){
  if (is.character(regi)){
    if (is.null(attr(x, 'regions'))){
      stop('No attribute with region names')
    } else {
      regi <- match(regi, attr(x, 'regions'))    
    } 
  } else if (is.logical(regi)){
    if (is.null(attr(x, 'regions'))){
      stop('No attribute with region names')
    } else if (length(attr(x, 'regions')) != length(regi)){
      stop('length of logical region selector does not match number of regions')
    } else {
      regi <- which(regi)
    }
  } 
  ndim <- length(dim(x))
  xtmp <- collapse2mat(aperm(x, c(ndim - 1, setdiff(1:ndim, ndim-1))), first=TRUE)
  outdims <- dim(x)[-(ndim-1)]
  if (length(regi) > 1){    
    outdims <- c(length(regi), outdims)
  }
  if (length(outdims) == 1) {
    out <- xtmp[regi,,drop=F]
  } else if (length(outdims) == ndim) {
    out <- aperm(array(xtmp[regi,], outdims), c(setdiff(1:ndim, c(1,ndim)),1, ndim))
  } else if (length(outdims) == ndim - 1) {
    out <- array(xtmp[regi,], outdims)
  }
  attnames <- setdiff(names(attributes(x)), c('dim', 'dimnames'))
  for (atn in attnames){
    attr(out, atn) <- attr(x, atn)
    if (atn %in% c('lat', 'lon')){
      if (length(attr(x, atn)) == dim(x)[ndim - 1]){
        attr(out, atn) <- attr(x, atn)[regi]
      }
    }
  }
  if (length(outdims) == 1){
    adimnames <- attr(x, 'dimnames')
    if (!is.null(adimnames[[length(adimnames) - 1]])) adimnames[[length(adimnames) - 1]] <- adimnames[[length(adimnames) - 1]][regi]
  } else {
    adimnames <- attr(x, 'dimnames')[-2]
  }
  attr(out, 'dimnames') <- adimnames
  if ('regions' %in% names(attributes(x))) attr(out, 'regions') <- attr(x, 'regions')[regi]
  invisible(out)
}
