#' Select specified years
#' 
#' Function to select a range of year from object of class 'NetCDF'.
#' 
#' @param x object of class 'NetCDF'
#' @param startyear first year to be selected
#' @param endyear last year to be selected
#' 
#' @keywords utilities
#' @export
select_time   <- function(x, startyear, endyear){
  x.time  <- attr(x, "time")
  if (!is.null(x.time)){
    if (min(x.time) < 100 && startyear > 100){
      # do not subset times for control runs
      out <- x
    } else {
      ndim    <- length(x.time)
      x.dims  <- dim(x)
      
      if (length(x.dims) < 2){
        x.dims <- length(x)
        
        # set up logical index
        tim.ind <- x.time >= startyear & x.time < (endyear + 1) & !is.na(x.time)
        if (sum(tim.ind) > 0){
          out <- x[tim.ind]
          attnames <- setdiff(names(attributes(x)), 'dim')
          for (atn in attnames){
            attr(out, atn) <- attr(x, atn)
          }
          attr(out, "time") <- x.time[tim.ind]
        }               
        
      } else {
        # find the time dimension
        time.i  <- which(x.dims == ndim)
        if ( length(time.i) != 1){
          warning("Dimensions probably messed up in select_time")
          time.i <- time.i[length(time.i)]
        }
        n.before<- prod(x.dims[1:(time.i-1)])
        
        # set up logical index
        tim.ind <- x.time >= startyear & x.time < (endyear + 1) & !is.na(x.time)
        
        if (sum(tim.ind) > 0){
          # output dimensions
          out.dim <- x.dims
          out.dim[time.i] <- sum(tim.ind)
          
          index   <- rep(tim.ind, each=n.before, length.out=length(x))
          
          out <- array(x[index], out.dim)
          attnames <- setdiff(names(attributes(x)), 'dim')
          for (atn in attnames){
            attr(out, atn) <- attr(x, atn)
          }
          attr(out, "time") <- x.time[tim.ind]
        }
      }
    }
  } else {
    out <- x
  }
  invisible(out)
}
