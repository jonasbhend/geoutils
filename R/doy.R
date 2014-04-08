#' Compute the day of the year from date
#'
#' @param x Rdates to compute the day of the year
#' @keywords utilities
#' @export
doy <- function(x) {
  out <- rep(NA, length(x))
  iix <- !(!is.finite(x) | is.na(x))
  if (any(iix)){
    year <- format(x[iix], '%Y')
    refdate <- as.Date(paste(year, '01', '01', sep='-'))
    out[iix] <- julian(x[iix]) - julian(refdate) + 1
  }
  return(out)
}
