#' Compute p-values when controlling the False Discovery Rate
#'
#' Function to compute adjusted critical p-values and the fraction of
#' test for which the alternative hypothesis might be true according
#' to Ventura et al. (2004).
#'
#' @param p p-values of multiple tests
#' @param q probability for which FDR should be controlled
#' @param log logical, are p-values log-transformed?
#' @param tight logical, should tighter FDR control be used
#' @keywords utilities
#' @references Ventura, V., C.J. Paciorek, and J.S. Risbey (2004). Controlling the proportion of falsely rejected hypotheses when conducting multiple tests with climatological data. Journal of Climate, 17 (22).
#' @export
fdrfun <- function(p, q=0.05, log=TRUE, tight=TRUE){
  if (log){
    np <- exp(p)
    bonferroni <- log(q / sum(!is.na(p)))
  } else {
    np <- p
    bonferroni <- q / sum(!is.na(p))
  }
  if (tight){
    ## compute tighter upper control
    xi <- seq(0.6, 0.9999, length=20)[-20]
    fp <- function(x) apply(as.matrix(x), 1, function(y) mean(np < y, na.rm=T))
    a <- 1/20 * sum(pmax(0, (fp(xi) - xi)/(1-xi)))
    qmod <- q/(1 - a)
  } else {
    qmod <- q
    a <- NA
  }
  
  if (log){
    iq <-log(qmod*seq(sort(p))/ sum(!is.na(p)))
  } else {
    iq <- qmod*seq(sort(p))/sum(!is.na(p))
  }
  
  if (any(sort(p) < iq)){
    pa <- max(bonferroni, sort(p)[max(which(sort(p) < iq))])
  } else {
    pa <- bonferroni
  }
  
  return(list(pcrit=pa, ahat=a))
}


