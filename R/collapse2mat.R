#' Converts array to matrix
#' 
#' Function to convert an array with at least 2 dimensions to a matrix.
#' Either the first dimension (as rows) or the last (as columns) will be
#' retained.
#'
#' @param x array with at more than two dimensions
#' @param first logical, should first or last dimension be retained?
#'
#' @keywords utilities
#' @export
#' @examples
#' ## should be 3, 20
#' dim(collapse2mat(array(NA, c(3,5,4)), first=TRUE))
collapse2mat <- function(x, first=FALSE){
  dims <- dim(x)
  if (length(dims) <= 1) dims <- c(1, length(x))
  if (first){
    xout <- array(x, c(dims[1], length(x)/dims[1]))
  } else {
    xout <- array(x, c(prod(dims[2:length(dims) - 1]), dims[length(dims)]))
  }
  xout
}
