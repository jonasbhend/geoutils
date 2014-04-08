#' New \code{open.ncdf} for CF compatibility
#'
#' This function replaces the \code{\link[ncdf]{open.ncdf}} in the \code{\link{ncdf}}
#' package to ensure compliance with CF standard with regard to the missing
#' value attribute.
#' 
#' @param con Name of the existing netCDF file to be opened.
#' @param write If \code{FALSE} (default), then the file is opened read-only.
#' If \code{TRUE}, then writing to the file is allowed.
#' @param readunlim When invoked, this function reads in the values of all
#' dimensions from the associated variables. This can be slow for a large file
#' with a long unlimited dimension. If set to \code{FALSE}, the values for the
#' unlimited dimension are not automatically read in (they can be read in 
#' later, manually, using \code{\link[ncdf]{get.var.ncdf}()}).
#' @param ... Arguments passed to or from other methods.
#' @param verbose	If \code{TRUE}, then messages are printed out during 
#' execution of this function.
#' 
#' @export
open.ncdf <- function( con, write=FALSE, readunlim=TRUE, verbose=FALSE, ... ) {
  
  if( verbose )
    print(paste("open.ncdf: entering, version=",version.ncdf()))
  
  rv <- list()
  
  if( write )
    rv$cmode <- 1
  else
    rv$cmode <- 0
  
  rv$id    <- -1
  rv$error <- -1
  rv <- .C("R_nc_open",
           as.character(con),
           as.integer(rv$cmode),
           id=as.integer(rv$id),
           error=as.integer(rv$error),
           PACKAGE="ncdf")
  if( rv$error != 0 ) 
    stop(paste("Error in open.ncdf trying to open file",con))
  if( verbose )
    print(paste("open.ncdf: back from call to R_nc_open, ncid=",rv$id))
  
  #-------------------------------------------------
  # Now we make our elaborate ncdf class object
  #-------------------------------------------------
  nc       <- list()
  nc$id    <- rv$id
  attr(nc,"class") <- "ncdf"
  
  #---------------------------------------
  # Get general information about the file
  #---------------------------------------
  if( verbose )
    print("open.ncdf: getting general info about the file")
  rv            <- list()
  rv$ndims      <- -1
  rv$nvars      <- -1
  rv$natts      <- -1
  rv$unlimdimid <- -1
  rv$error      <- -1
  rv <- .C("R_nc_inq",
           as.integer(nc$id),
           ndims=as.integer(rv$ndims),
           nvars=as.integer(rv$nvars),
           natts=as.integer(rv$natts),
           unlimdimid=as.integer(rv$unlimdimid),
           error=as.integer(rv$error),
           PACKAGE="ncdf")
  if( rv$error != 0 ) 
    stop(paste("R_nc_inq returned error on file",con,"!"))
  if( verbose )
    print(paste("open.ncdf: back from call to R_nc_inq for id=",nc$id,"   ndims=",rv$ndims,
                "   nvars=",rv$nvars,"  natts=",rv$natts,"  umlimdimid=",rv$unlimdimid))
  nc$ndims        <- rv$ndims
  nc$natts        <- rv$natts
  nc$unlimdimid   <- rv$unlimdimid + 1	# Change from C to R convention
  nc$filename     <- con
  nc$varid2Rindex <- rep(0,rv$nvars)	
  nc$writable     <- write
  
  #-----------------------------------------------
  # Get all the dimensions that this file has.
  # Get their values as well (for caching).  
  #-----------------------------------------------
  nc$dim   <- list()
  dimnames <- character()
  for( i in 1:nc$ndims ) {	
    if( verbose )
      print(paste("open.ncdf: getting dim info for dim number",i))
    #-----------------------------------------------------------
    # As a general note, the function dim.inq.ncdf does NOT
    # return a full-fledged "dim.ncdf" object. It returns 
    # only a subset of fields that directly correspond to 
    # a low-level, netCDF dimension in a file.  We now fill
    # out the rest of the fields to make it into a real dim.ncdf 
    # object.
    #-----------------------------------------------------------
    d          <- dim.inq.ncdf(nc,i)
    d$id	   <- i
    d$dimvarid <- varid.inq.ncdf(nc,d$name)
    if( verbose )
      print(paste(".....dim name is",d$name,"  len=",d$len,"     dimvarid=",d$dimvarid))
    if( d$dimvarid == -1 ) {	# No dimvar for this dim
      d$vals  <- 1:d$len
      d$units <- ""
      d$create_dimvar <- FALSE	# in case this dim is passed to create.ncdf()
    }
    else {	
      # This dim has a dimvar -- get its properties
      if( verbose )
        print(paste("open.ncdf: getting dimvar info for dim ",d$name))
      attv <- att.get.ncdf( nc, d$dimvarid, "units" )
      if( attv$hasatt )
        d$units <- attv$value
      else
        d$units <- ""
      if( d$unlim && (! readunlim)) # Is unlimited, don't read vals, too slow
        d$vals <- rep(NA,d$len)
      else			# Otherwise, read vals
        d$vals <- get.var.ncdf( nc, forcevarid=d$dimvarid, verbose=verbose )
      d$create_dimvar <- TRUE		# in case this dim is passed to create.ncdf()
      if( verbose )
      {
        print("------------------------------")
        print("Here is new dim:")
        print(paste("name=",d$name,"  len=",d$len,"   unlim=",d$unlim,"   id=",d$id,"   dimvarid=",d$dimvarid,"   units=",d$units))
        print("------------------------------")
      }
    }
    attr(d,"class") <- "dim.ncdf"	# Is a complete dim.ncdf object now
    nc$dim[[i]] <- d
    dimnames[i] <- d$name
    if( verbose )
      print(paste(".......open.ncdf: done processing dim ",d$name))
  }
  attr(nc$dim,"names") <- dimnames
  if(verbose) {
    print("open.ncdf: setting dim$<names> to:")
    print(dimnames)
  }
  
  #-------------------------------------------
  # Get all the vars that this file has.  Note
  # that dimvars are NOT included in the count
  # of vars!!
  #-------------------------------------------
  if( verbose )
    print(paste("open.ncdf: getting var info.  # vars (INCLUDING dimvars)=",rv$nvars))
  nc$nvars <- 0
  nc$var   <- list()
  varnames <- character()
  for( i in 1:rv$nvars ) {
    name <- varname.inq.ncdf( nc, i )
    if( dimid.inq.ncdf( nc, name ) == -1 ) {	# Only process if NOT a dimvar
      if( verbose )
        print(paste("open.ncdf: will process varid=",i,"  name=",name))
      #--------------------------------------
      # No dim with same name as this var, so
      # this var must NOT be a dimvar.
      #--------------------------------------
      v <- var.inq.ncdf( nc, i )
      attr(v,"class") <- "var.ncdf"
      nc$nvars <- nc$nvars + 1
      #--------------------
      # Get this var's dims
      #--------------------
      v$dims   <- list()
      varunlim <- FALSE
      if( v$ndims > 0 ) {
        for( j in 1:v$ndims ) {
          v$dim[[j]] <- nc$dim[[ v$dimids[j] ]]
          if( v$dim[[j]]$unlim )
            varunlim <- TRUE
          v$varsize <- append(v$varsize, v$dim[[j]]$len)
        }
      }
      v$unlim <- varunlim
      
      #----------------------------------------
      # Get this var's missing value, or set to
      # a default value if it does not have one
      #----------------------------------------
      mv <- att.get.ncdf( nc, i, "_FillValue" )
      if( mv$hasatt ) {
        v$missval <- mv$value
      }
      else {
        mv <- att.get.ncdf( nc, i, "missing_value" )
        if( mv$hasatt )
          v$missval <- mv$value
        else if( (v$prec=="float") || (v$prec=="double"))
          v$missval <- default.missval.ncdf()
        else
          v$missval <- NA
      }
      
      #-------------------------------------------
      # Get add_offset and scale_factor attributes 
      #-------------------------------------------
      ao <- att.get.ncdf( nc, i, "add_offset" )
      if( ao$hasatt ) {
        v$hasAddOffset <- TRUE
        v$addOffset    <- ao$value
      }
      else
        v$hasAddOffset <- FALSE
      sf <- att.get.ncdf( nc, i, "scale_factor" )
      if( sf$hasatt ) {
        v$hasScaleFact <- TRUE
        v$scaleFact    <- sf$value
      }
      else
        v$hasScaleFact <- FALSE
      
      nc$var[[nc$nvars]] <- v
      nc$varid2Rindex[i] <- nc$nvars 
      varnames <- append(varnames,v$name)
      if( verbose ) {
        print("-----------------------")
        print("Here is new var:")
        print(paste("name=",v$name,"  id=",v$id,"   ndims=",v$ndims,"   prec=",v$prec))
        print("size=")
        print(v$size)
        print("dimids=")
        print(v$dimids)
      }
    }
  }
  attr(nc$var,"names") <- varnames
  
  if( verbose )
    print(paste("open.ncdf: leaving for ncid=",nc$id))
  
  return(nc)
}
