#' Select a region
#' 
#' Function to select regions (2nd last dimension in NetCDF array).
#' 
#' @param x input array of class NetCDF
#' @param regi index of region to select
#' 
#' @keywords utilities
#' @export
select_region <- function(x, regi=1){
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
    adimnames[[length(adimnames) - 1]] <- adimnames[[length(adimnames) - 1]][regi]
  } else {
    adimnames <- attr(x, 'dimnames')[-2]
  }
  attr(out, 'dimnames') <- adimnames
  if ('regions' %in% names(attributes(x))) attr(out, 'regions') <- attr(x, 'regions')[regi]
  invisible(out)
}
