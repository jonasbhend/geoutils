#' Compute distance between points on the globe
#' 
#' This function computes distances between an arbitrary number of points on the globe
#'
#' @param lon,lat points from which distance is measured (in degree)
#' @param lon0,lat0 points to which distance is measured (in degree)
#' @param Rearth Earth's radius (defaults to 6371km)
#' @return Distance between points in km
#' @keywords utilities
#' @export
compute_dist <- function(lon, lat, lon0, lat0, Rearth=6371){
  lo <- lon/180*pi
  la <- lat/180*pi
  lo0 <- lon0/180*pi
  la0 <- lat0/180*pi
  tmp <- outer(sin(la), sin(la0), '*') + outer(cos(la), cos(la0), '*')*outer(lo, lo0, function(x,y) cos(x - y))
  # deal with rounding errors
  tmp[abs(tmp) >= 1] <- round(tmp[abs(tmp) >= 1])
  if (ncol(tmp) == 1) tmp <- as.vector(tmp)
  out <- Rearth*acos(tmp)
  out
}
