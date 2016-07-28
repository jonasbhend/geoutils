#' Read in NetCDF data
#'
#' This function reads in NetCDF data and stores it in an object with
#' class 'NetCDF'
#'
#' @param filename path to input NetCDF file
#' @param varname variable name to read in (one at a time)
#' @param mask vector (logical) with lon-lat indices to be read in
#' @param mulc Multiplicative factor to be multiplied with NetCDF data (see details)
#' @details
#' The multiplicative factor can be of length 1 or any subset of the time
#' dimension of the input data so that seasonal data can be converted from
#' units per month to units per season with different multiplicative constants
#' per season.
#' 
#' Units attribute will not get converted if multiplicative constant is different
#' from unity.
#' 
#' @keywords utilities
#' @examples
#' tas <- readNetCDF(system.file("extdata", "annual_CRUTEMv3_1961-90.nc", package="geoutils"), varname="temp")
#' names(attributes(tas))
#' dim(tas)
#' @export
#' 
#' 
readNetCDF <- function(filename, varname=NULL, mask=NULL, mulc=1){
  nc <- open.ncdf(filename)
  on.exit(close.ncdf(nc))
  if (is.null(varname)) {
    varnames <- names(nc$var)
    if (length(grep('tim', varnames)) > 0){
      varnames <- varnames[-grep('tim', varnames)]
    }
    varname <- varnames[1]
  }
  dimn <- sapply(nc$var[[varname]]$dim, function(x) x$name)
  with.time <- FALSE
  timeid <- grep('tim', dimn)
  if (length(timeid) > 0){
    if (nc$dim[[dimn[timeid]]]$len > 1){
      with.time <- TRUE
    }
    dimid <- nc$var[[varname]]$dimid[-timeid]
  } else {
    dimid <- nc$var[[varname]]$dimid
  }
  dims <- list()
  if (length(dimid) > 1){
    dims[[names(nc$dim)[dimid[1]]]] <- rep(nc$dim[[dimid[1]]]$vals,
                                           length(nc$dim[[dimid[2]]]$vals))
    dims[[names(nc$dim)[dimid[2]]]] <- rep(nc$dim[[dimid[2]]]$vals,
                                           each=length(nc$dim[[dimid[1]]]$vals))
  } else {
    dims[[names(nc$dim)[dimid[1]]]] <- nc$dim[[dimid[1]]]$vals
    ## check whether longitudes and latitudes are available
    if (all(c('lon', 'lat') %in% names(nc$var))){
      dims[['lon']] <- get.var.ncdf(nc, 'lon')
      dims[['lat']] <- get.var.ncdf(nc, 'lat')
    }
  }
  if (with.time){
    tim <- ncdf_times(nc, timename=dimn[timeid])
    ind <- rep(TRUE, length(tim))
    nseas <- 12/median(table(as.numeric(format(tim, '%Y'))))
    if (nseas == 3){
      if (length(tim) %% (12/nseas) == 1){
        if (as.numeric(format(tim, '%m'))[1] == 11){
          ind[1] <- FALSE
        } else {
          ind[length(ind)] <- FALSE
        }
      }
      if (all(as.numeric(format(tim, '%d')) > 25)){
        dims$time <- (as.numeric(format(tim, '%Y')) + (as.numeric(format(tim, '%m')) - 1)/12)[ind]
      } else {
        dims$time <- (as.numeric(format(tim, '%Y')) + (as.numeric(format(tim, '%m')) - 2)/12)[ind]
      }
      ## set the first value to missing (consisting of two months only)
      ind <- which(ind)
      if (length(tim) > 4) ind[1] <- NA
    } else if (nseas == 12){
      dims$time <- as.numeric(format(tim, '%Y'))
    } else {
      dims$time <- as.numeric(format(tim, '%Y')) + (as.numeric(format(tim, '%m')) - 1)/12
    }
    ## read in the data
    out <- collapse2mat(get.var.ncdf(nc, varname), first=FALSE)
    ## adjust scaling (works with different values per season)
    ## make sure that mulc has correct length
    mulc <- rep(mulc, length=12 / nseas)
    out <- out * rep(mulc, each=nrow(out), length=length(out))
    if (!is.null(mask) & length(mask) == nrow(out)) {
      out <- out[mask,ind, drop=F]
      if (!is.null(mask)) {
        for (i in which(sapply(dims, length) == length(mask) & names(dims) != 'time')){
          dims[[i]] <- dims[[i]][mask]
        }
      }
    } else {
      out <- out[,ind, drop=F]
    }
  } else {
    ## read in the data
    out <- as.vector(get.var.ncdf(nc, varname))*mulc[1]
    if (!is.null(mask) & length(mask) == length(out)) {
      out <- out[mask]
      for (i in which(sapply(dims, length) == length(mask))){
        dims[[i]] <- dims[[i]][mask]
      }
    }
    out <- as.matrix(out)
    ## fix to get the time attribute back in for single record files
    if (length(timeid) > 0){
      tim <- ncdf_times(nc, timename=dimn[timeid])
      if (all(as.numeric(format(tim, '%d')) > 25)){
        dims$time <- (as.numeric(format(tim, '%Y')) + (as.numeric(format(tim, '%m')) - 1)/12)
      } else if (format(tim, '%m') == '01' & format(tim, '%d') == '01') {
        dims$time <- as.numeric(format(tim, '%Y'))
      } else {
        dims$time <- (as.numeric(format(tim, '%Y')) + (as.numeric(format(tim, '%m')) - 2)/12)
      }
    }
  }
  ## time check for kaputt times in HadGEM2-CC
  if (any(diff(dims$time) > 100*min(diff(dims$time)))) dims$time <- seq(dims$time[1], length=length(dims$time), by=diff(dims$time)[1])
  ## check for kaputt values
  
  for (n in names(dims)){
    attr(out, n) <- dims[[n]]
  }
  if (att.get.ncdf(nc, 0, 'forcing')$hasatt){
    attr(out, 'forcing') <- att.get.ncdf(nc, 0, 'forcing')$value
  }
  attr(out, 'units') <- nc$var[[varname]]$units
  if (any(names(nc$var) %in% c('region', 'regions', 'cluster', 'clusters'))){
    nname <- names(nc$var)[names(nc$var) %in% c('region', 'regions', 'cluster', 'clusters')]
    attr(out, 'regions') <- get.var.ncdf(nc, nname[1])
  } else if (att.get.ncdf(nc, 0, 'region_names')$hasatt) {
    attr(out, 'regions') <- strsplit(att.get.ncdf(nc, 0, 'region_names')$value, ', ')[[1]]
  }
  neffatt <- att.get.ncdf(nc, 0, 'number_of_effective_grid_boxes')
  if (neffatt$hasatt) attr(out, 'neff') <- as.numeric(strsplit(neffatt$val, ',')[[1]][-1])
  class(out) <- 'NetCDF'
  return(out)
}
