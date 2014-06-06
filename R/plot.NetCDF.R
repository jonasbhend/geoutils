#' Function to plot objects of class 'NetCDF'
#' 
#' This function allows plotting line plots and maps of NetCDF objects using
#' the respective attributes.
#' 
#' @param x input object of class NetCDF (max. 2-d)
#' @param type of plot ('ts' for lines, 'mean' and 'trend' for maps)
#' @param ti index of time step to be plotted (map only)
#' @param si index of spatial grid to be plotted
#' @param levs levels to be used for contouring
#' @param col colours to be used for contouring or lines
#' @param na.col colour to be used for missing values (spatial plot only)
#' @param pt.cex expansion factor for points
#' @param symmetric logical, should contouring be symmetric about zero
#' @param cut logical, should contouring be cut to range of data
#' @param add logical, should plot be added
#' @param xlab,ylab axes labels
#' @param xlim,ylim range of axes
#' @param xaxs,yaxs type of axes
#' @param lty line type for line plots
#' @param lwd line width for line plots
#' @param nlevs number of contour interval (input to \code{pretty})
#' @param colramp colour ramp to be used with colourramp function (default)
#' @param seas logical, should seasons be plotted individually?
#' @param ... additional arguments passed to \code{plot.default}
#' 
#' @keywords plot
#' @examples
#' tas <- readNetCDF(system.file("extdata", "annual_CRUTEMv3_1961-90.nc", package="geoutils"), varname="temp")
#' ## average over all time steps
#' plot(tas)
#' map2()
#' 
#' ## also plot the trends
#' plot(tas, type='trend', sym=TRUE)
#' map2()
#' 
#' ## plot the same with legend
#' layout(matrix(1:2, 2, 1), height=c(5, lcm(2.5)))
#' par(mar=c(2,2,1,1))
#' tmp <- plot(tas, type='trend', sym=TRUE, na.col='lightgrey')
#' ## if pt.cex is not set, the default pt.cex will be displayed
#' map2()
#' plot_colourbar(tmp, units='K/a')
#' 
#' 
#' @export
plot.NetCDF <- function(x, type=if (nrow(x) == 1) 'ts' else 'mean', ti=1:dim(x)[length(dim(x))], si=1, levs=NULL, col=NULL, na.col=NULL, pt.cex=NULL, symmetric=FALSE, cut=TRUE, add=FALSE, xlab='', ylab='', xlim=NULL, ylim=NULL, xaxs='i', yaxs='i', lty=seq(along=si), lwd=3, nlevs=12, colramp='redblue', seas=F, ...){
  ###### time series plots #############################
  if (type == 'ts'){
    xtmp <- select_region(x, si)
    tim <- if ('time' %in% names(attributes(x))) attr(x, 'time') else 1:ncol(x)
    ntim <- length(tim)
    xtmp <- collapse2mat(xtmp, first=F)
    nseas <- median(table(floor(tim)))
    if (is.null(col)){
      cols <- seq(1,nseas)
    } else {
      cols <- rep(col, length=nrow(xtmp)*nseas)
    }
    if (!add) plot(tim, xtmp[1,], type='n', xlab=xlab, ylab=ylab, xlim=xlim, ylim=if (is.null(ylim)) range(xtmp, na.rm=T) else ylim, ...)
    if (seas){
      for (i in 1:nseas){
        for (j in seq(along=si)){
          ind <- seq(i,ntim,nseas)
          lines(tim[ind], xtmp[j,ind], lwd=lwd, col=cols[i], lty=rep(lty, length=length(si))[j])        
        }
      }      
    } else {
      for (j in seq(along=si)){
        lines(tim, xtmp[j,], lwd=lwd, col=cols[1], lty=rep(lty, length=length(si))[j])
      }
    }
    ###### spatial plots #############################
  } else {
    if ('lon' %in% names(attributes(x)) & 'lat' %in% names(attributes(x))){
      lon <- attr(x, 'lon')
      lat <- attr(x, 'lat')
    } else {
      ## get all attributes
      allatt <- lapply(names(attributes(x)), function(y) attr(x, y))
      iatt <- which(sapply(allatt, length) == nrow(x) & sapply(allatt, is.numeric))
      if (length(iatt) >= 2){
        lon <- allatt[[iatt[1]]]
        lat <- allatt[[iatt[2]]]
      } else{
        stop('No spatial coordinates found')
      }
    } 
    xtmp <- x
    if (dim(xtmp)[1] == 1 & length(dim(xtmp)) == 3) xtmp <- xtmp[1,,]
    if (length(dim(xtmp)) > 2){
      xtmp <- apply(xtmp, length(dim(x)) - 1:0, mean, na.rm=T)
    } 
    if (type == 'mean'){
      xout <- apply(xtmp[,ti, drop=F], 1, mean, na.rm=T)
    } else if (type == 'trend' & length(ti) > 1) {
      xtim <- attr(x, 'time')
      xout <- apply(xtmp, 1, function(x) if (any(!is.na(x))) lm(x[ti] ~ xtim[ti])$coef[2] else NA)
    }
    if (is.null(levs)) {
      if (symmetric & cut){
        levs <- pretty(xout, nlevs)
        levs2 <- sort(unique(c(levs, - levs)))
        col <- colourramp(length(levs2) - 1, ramp=colramp)
        col <- col[levs2[-1] %in% levs[-1]]  
        rm(levs2)    
      } else if (symmetric & ! cut) {
        levs <- pretty(c(xout, -xout), nlevs)
      } else {
        levs <- pretty(xout, nlevs)
      }
    }
    if (is.null(col)) col <- colourramp(length(levs) -1, ramp=colramp)
    if (is.null(xlim)){
      xlim <- range(lon)+ c(-0.5,0.5)*median(abs(diff(lon)[diff(lon) != 0]))
    }
    if (is.null(ylim)) {
      ylim <- range(lat) + c(-0.5, 0.5)*median(abs(diff(lat)[diff(lat) != 0]))
    }
    if (min(lon) < xlim[1]){
      lon[lon < xlim[1]] <- lon[lon < xlim[1]] + 360
    }
    if (max(lon) > xlim[2]){
      lon[lon > xlim[2]] <- lon[lon > xlim[2]] - 360
    }
    if (!add) plot(lon, lat, type='n', xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim, xaxs=xaxs, yaxs=yaxs, ...)
    ## figure out size of plot symbol if not set
    if (is.null(pt.cex)){
      ## get width and height
      wh <- diff(par('usr'))[c(1,3)]
      dimperpoint <- wh / c(length(unique(lon)), length(unique(lat)))
      ## standard size of character
      psize <- par('cxy')
      ## point size
      pt.cex <- max(2*dimperpoint / psize)
      print(paste('pt.cex =', round(pt.cex,2)))
    }    
    points(lon,lat,col=col[cut(xout, levs, include.lowest=TRUE)], pch=15, cex=pt.cex)
    ## add in colour for missing values
    if (!is.null(na.col) & any(is.na(xout))) points(lon[is.na(xout)], lat[is.na(xout)], pch=15, cex=pt.cex, col=na.col)
    out <- list(lev=levs, col=col, sea.col=na.col)
    invisible(out)
  }
}
