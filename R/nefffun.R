#' Compute effective degrees of freedom for autocorrelated time series
#'
#' @param N number of elements in input
#' @param acf autocorrelation coefficients (up to lag 2, higher lags not implemented)
#' @keywords utilities
#' @export
nefffun <- function(N, acf){
  if (length(acf) < N - 1){
    if (length(acf) == 1){
      acf <- acf[1] ** seq(1,N-1)
    } else if (length(acf) == 2) {
      phi1 <- acf[1] * (1 - acf[2])/(1 - acf[1]**2)
      phi2 <- (acf[2] - acf[1]**2)/(1 - acf[1]**2)
      for (i in 3:N) acf <- c(acf, phi1*acf[i-1] + phi2 * acf[i-2])
    } else {
      stop('Higher lags not implemented in nefffun')
    }
  }
  neff <- N/(1 + 2*sum((1-seq(1,N-1)/N)*acf[seq(1,N-1)]))
  #if (neff > N) {
  #  warning('Reset neff to N')
  #  neff <- N
  #}
  return(neff)
}
