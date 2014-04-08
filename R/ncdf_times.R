#' Get time attribute from NetCDF file and convert it to years (or Rdate)
#'
#' This function reads the NetCDF time attribute and converts it to years
#' with decimal fraction for monthly data. For sub-monthly data the conversion
#' is to the Rdate format (only works with standard calendar!).
#'
#' @param nc object with link to NetCDF file (from \code{\link{open.ncdf}})
#' @param timename dimension name with time (optional)
#' @param as.Rdate logical, should output be converted to Rdate?
#' @param force logical, force Rdate conversion for monthly times?
#' @param tz time zone for time zone support (experimental)
#' @keywords utilities
#' @export
#' 
ncdf_times <- function(nc, timename=NULL, as.Rdate=TRUE, force=TRUE, tz="CET") {
  ## this function converts netcdf times to the 
  ## R date-time format or to the udunits dates
  ## you can choose to switch to uduints format
  ## for gregorian calendars by setting as.Rdate
  ## to FALSE
  ## non-gregorian calendar dates are output using 
  ## udunits date format
  
  ## you can force to get back an R date format, even
  ## if the calendar used is not gregorian using
  ## force=T (may return udunits dates if conversion
  ## is not successful)
  
  ## WARNING: time zones are not fully supported yet
  
  ## check whether udunits is available
  .udunitsInclude     <- FALSE
  if (any(.packages() == "udunits") & class(try(utInit(), silent=T)) != "try-error"){
    .udunitsInclude <- TRUE
  }
  if (is.null(timename)){
    timei   <- which(names(nc$dim) %in% c("time", "TIME", "tim", "TIM"))
  } else {
    timei <- which(names(nc$dim) == timename)
  }
  units   <- nc$dim[[timei]]$units
  refdate <- strsplit(units, " ")[[1]]
  refdate <- refdate[grep('-', refdate)]
  ## debug reference date
  if (substr(refdate, nchar(refdate) - 1, nchar(refdate)) == '00') {
    rtmp <- strsplit(refdate, '-')[[1]]
    refdate <- paste(c(rtmp[-length(rtmp)], '01'), collapse='-')
    rm(rtmp)
  }
  vals    <- nc$dim[[timei]]$vals
  tmp     <- att.get.ncdf(nc, names(nc$dim)[timei], "calendar")
  if (tmp$hasatt) {
    calendar <- tmp$value
  } else {
    calendar <- "standard"
    ## print(paste("Warning: Calendar is missing in", nc$filename))
  }    
  if (calendar == "proleptic_gregorian" || calendar == "gregorian") calendar <- "standard"
  
  if (as.Rdate){
    
    if (charmatch("hours since", units, nomatch=0) | 
          charmatch("minutes since", units, nomatch=0) |
          charmatch("seconds since", units, nomatch=0) & calendar == 'standard') {
      
      mul <- 1
      ref.txt     <- substr(units, 15,33)
      if (charmatch("minutes", units, nomatch=0)) mul     <- 60
      if (charmatch("hours", units, nomatch=0)) {
        mul     <- 3600
        ref.txt <- substr(units, 13,31)
      }
      
      times       <- vals * mul
      if (nchar(ref.txt) == 19){
        ref         <- as.POSIXct(ref.txt, tz)
      } else {
        ref         <- as.POSIXct(paste(ref.txt, "00", sep=":"), tz)
      }   
      time        <- as.Date(ref + times)
      
    } else if (charmatch("days since", units, nomatch=0) & calendar == 'standard'){
      time        <- as.Date(refdate, "%Y-%m-%d") + vals
    } else if (charmatch("days since", units, nomatch=0) &
                 calendar %in% c('365_day', 'noleap', '360_day')) {
      if (calendar == '365_day' || calendar == 'noleap'){
        vals <- round(vals/365*365.24*2)/2
        time <- as.Date(refdate, "%Y-%m-%d") + vals         
      } else if (calendar == '360_day'){
        vals <- round(vals/360*365.24*2)/2
        time <- as.Date(refdate, "%Y-%m-%d") + vals        
      } 
    } else if (charmatch("months since", units, nomatch=0)) {
      ref.yr      <- as.numeric(format(as.Date(refdate), '%Y'))
      ref.month   <- as.numeric(format(as.Date(refdate), '%m'))
      ref.day     <- as.numeric(format(as.Date(refdate), '%d'))
      if (is.null(ref.day) | ref.day == 0) ref.day <- 1
      
      month       <- floor((vals+ref.yr*12 + ref.month-1) %% 12) + 1
      year        <- floor((vals+ref.yr*12 + ref.month-1)/12) 
      
      time        <- as.Date(ISOdate(year, month, ref.day))
    } else if (charmatch("years since", units, nomatch=0)) {
      unit.tmp    <- paste(strsplit(units, " ")[[1]][3:4], collapse=" ")
      ## ref.yr      <- substr(units, 13,16)
      ## ref.month   <- as.numeric(substr(units, 18,19))
      ref.yr      <- as.numeric(format(as.Date(unit.tmp), "%Y"))
      ref.month   <- as.numeric(format(as.Date(unit.tmp), "%m"))
      if (is.null(ref.month)) ref.month <- 1
      ##ref.day     <- as.numeric(substr(units, 21,22))
      ref.day     <- as.numeric(format(as.Date(unit.tmp), "%d"))
      if (is.null(ref.day)) ref.day <- 1
      
      year        <- floor(vals)
      month       <- floor((vals*12)%%12)
      
      time        <- as.Date(ISOdate(ref.yr + year, ref.month + month, ref.day))
    }  else if (charmatch("day as", units, nomatch=0)) {
      date        <- floor(vals)
      day         <- as.numeric(substr(date, nchar(date)-1, nchar(date)))
      if (all(day > 28)) date <- as.character(as.numeric(date) - max(day, na.rm=T) + 28)
      date        <- paste("000",date, sep="")
      date        <- substr(date, nchar(date)-7, nchar(date))
      time        <- as.Date(date, "%Y%m%d")
    } else {
      stop(paste("Can't deal with calendar", calendar))
    }
    
  } else {
    if (.udunitsInclude){
      time <- utCalendar(vals, units, calendar=calendar, style="array")  
      if (force){
        tmp  <- try(ISOdatetime(time$year, time$month, time$day, time$hour, 
                                time$minute, time$second, tz), silent=T)
        if (class(tmp)[1] != "try-error") time <- tmp
      }
    } else {
      warning("Package udunits cannot be loaded or initialized via utInit()")
    }
  }
  
  
  return(time)
  
}
