#' @name string
#' @aliases simplecap
#' @aliases firstcap
#' 
#' @title
#' Functions to modify text.
#'
#' @description
#' These functions modify text to capitalise the first letter of each
#' word (\code{simplecap}) or the very first letter (\code{firstcap}).
#'
#' @param x text, vector of type character
#' @keywords manip
#' @export
simplecap <- function(x) {
  sapply(strsplit(x, " "), function(s) paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" "))
}

#' @rdname string
#' @export
firstcap <- function(x) {
  xout <- paste0(toupper(substring(x, 1, 1)), substring(x, 2, pmax(2,nchar(x))))
  attributes(xout) <- attributes(x)
  return(xout)
}

