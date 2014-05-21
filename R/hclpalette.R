#' Plots a palette or two hcl hues (through white)
#' 
#' This function returns a number of colours that form an hcl palette
#' based on the two hues provided
#' 
#' @param x number of colours returned
#' @param hues two hues for the colour endpoints (in the interval from 0 to 360)
#' 
#' @keywords plot
#' @export
hclpalette <- function(x, hues=c(0,240)){
  hind <- c(rep(hues[1], ceiling(x/2)), rep(hues[2], floor(x/2)))
  ind <- seq(-4.5,4.5,length=x)
  ind2 <- seq(-1.5, 1.5, length=x)
  hcl(hind, l=99/5.5*(5.5 - abs(ind)), c=(1 - abs(1 - abs(ind2)))*70)
}
