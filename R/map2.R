#' Addition to plot proper map without borders
#' 
#' Plots map without borders and adds missing lakes and inland seas
#' 
#' @param interior logical, should country borders be plotted? defaults to FALSE
#' @param add logical, should map be added to existing plot? devaults to TRUE
#' @param ... Additional arguments passed to map
#' 
#' @example
#' map2(add=F)
#' 
#' @keywords utilities
#' @export
#' 
map2 <- function(interior=F, add=T, ...){
  map(interior=interior, add=add, ...)
  if (!interior){
    map(region='Caspian', add=TRUE, ...)
    map(region='.*lake', exact=F, add=T, ...)
  }
}
