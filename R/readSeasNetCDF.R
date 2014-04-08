#' Read in NetCDF data with NCI seasonal structure (one variable per season)
#' and collate to continuous seasonal time series
#' 
#' This function is a wrapper for the \code{\link{readNetCDF}} function to
#' read in data in the NCI seasonal format and convert it to the a continous
#' time series.
#'
#' @param filename path to input NetCDF file
#' @param varname variable name to read in (one at a time)
#' @param seas string with season names to be read in (see details)
#' @param mask vector (logical) with lon-lat indices to be read in
#' @param mulc Multiplicative factor to be multiplied with NetCDF data (see details)
#' @details
#' The variable names are collated from the \code{varname} and the season 
#' names (e.g. tas_djf). Please note that the sequence of seasons is 
#' unrestricted, and a time series with for example annual, summer and winter
#' values can be generated.
#' @keywords utilities
#' @export

readSeasNetCDF <- function(filename, varname=NULL, seas=c('djf', 'mam', 'jja', 'son'), mask=NULL, mulc=1){
  if (all(seas == '')){
    out <- readNetCDF(filename, varname=varname, mask=mask, mulc=mulc)
  } else {
    
    mulc <- rep(mulc, length=length(seas))
    
    ## get the position of a substring
    cpos <- function (str, sub, start = 1){
      lstr <- nchar(str)
      lsub1 <- nchar(sub) - 1
      if (start + lsub1 > lstr) 
        return(NA)
      else {
        str <- substring(str, start, lstr)
        str <- substring(str, 1:(lstr - lsub1), (1 + lsub1):lstr)
        p <- ((start:lstr)[str == sub])[1]
        if (is.na(p > 0)) 
          return(NA)
        else return(p)
      }
    }
    
    ## read in seasonal time series 
    interim <- list()
    ## get season names for exception handling
    nc <- open.ncdf(filename)
    snames <- names(nc$var)
    close.ncdf(nc)
    snames <- gsub(paste(varname, '_', sep=''), '', snames[grep(varname, snames)])
    for (si in seq(seas)){
      se <- seas[si]
      ## add a whole lot of exceptions for e.g. Vietnam
      if (se %in% snames){
        interim[[si]] <- readNetCDF(filename, varname=paste(varname, se, sep='_'), mask=mask, mulc=mulc[si]) 
      } else {
        
        monstring <- 'jfmamjjasondjfmam'
        montxt <- c('january', 'february', 'march', 'april', 'may', 'june', 'july', 'august', 'september', 'october', 'november', 'december', 'january', 'february', 'march', 'april', 'may')
        dlen <- c(31,28,31,30,31,30,31,31,30,31,30,31, 90, 92, 92, 91, 181, 184)
        names(dlen) <- c(montxt[1:12],'djf', 'mam', 'jja', 'son', 'ndjfma', 'mjjaso')
        ss <- list()
        ss[[1]] <- montxt[cpos(monstring, se)]
        for (i in 2:nchar(se)) ss[[i]] <- montxt[cpos(monstring, se) + i - 1]
        ii <- lapply(ss, function(x) readNetCDF(filename, varname=paste(varname, x, sep='_'), mask=mask)*dlen[x])
        dtot <- sum(dlen[unlist(ss)])
        ## mask out first year if season crosses year
        if (any(unlist(ss) == 'january')){
          jani <- which(unlist(ss) == 'january')
          if (jani > 1){
            ii[seq(1,jani-1)] <- lapply(ii[seq(1,jani-1)], function(x) {
              xout <- x[,c(NA, seq(1, ncol(x) - 1))]
              attributes(xout) <- attributes(x)
              return(xout)})
          }
        }
        ## add all the months together
        interim[[si]] <- ii[[1]]
        for (i in seq(2, length(ii))) interim[[si]] <- interim[[si]] + ii[[i]]
        interim[[si]] <- interim[[si]] / dtot * mulc[si]
      }
    }
    ## convert to old season format
    out <- array(NA, dim(interim[[1]])*c(1,length(seas)))
    for (si in seq(seas)){
      out[,seq(si, ncol(out), length(seas))] <- interim[[si]]
    }
    ## add attributes
    attributes(out) <- c(attributes(out), attributes(interim[[1]])[! names(attributes(interim[[1]])) %in% c('dim', 'dimnames')])
    if (!is.null(rownames(interim[[1]]))) rownames(out) <- rownames(interim[[1]])
    attr(out, 'time') <- sort(as.vector(sapply(interim, attr, 'time'))) + (seq(seas) - 1)/length(seas)
    if (length(attr(out, 'time')) != ncol(out)) warning(paste('Time dimension messed up'))
  }
  return(out)
}
