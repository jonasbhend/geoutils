#' @name inset
#' @aliases reset_pos
#' 
#' @title functions to set position and size for insets
#' 
#' @description
#' \code{set_position} sets up a new plot at the specified position (in user
#' coordinates), \code{reset_position} resets the plot to the original setup.
#' 
#' @param inset.x,inset.y user coordinates for inset
#' @param fun Function to be evaluated for inset (e.g. \code{plot})
#' @param inset.size relative width and height of inset (compared to main plot)
#' @param inset.bg colour of background for inset
#' @param inset.mfcol logical, is original layout with multiple panels?
#' @param ... additional arguments passed to \code{fun}
#' 
#' @keywords utilities
#' @export
inset <- function(inset.x, inset.y, fun, inset.size=0.1, inset.bg=NULL, inset.mfcol=TRUE, ...){
    # sets the display for a plot centered at x,y
    oldpar    <- par(no.readonly=TRUE)
    on.exit(reset_pos(oldpar=oldpar, mfcol=inset.mfcol))
    lims        <- par("usr") 
    flims       <- par("plt")

    # expand size if necessary
    if (length(inset.size) %in% c(1,2)){
        inset.size    <- c(-0.5,0.5,-0.5,0.5)*rep(inset.size, each=2)
    } 
    
    if ( inset.x > lims[1] & inset.x < lims[2] & inset.y > lims[3] & inset.y < lims[4]){
        # change to plot coordinates
        figx    <- (inset.x-lims[1])/diff(lims[1:2])*diff(flims[1:2]) + flims[1]
        figy    <- (inset.y-lims[3])/diff(lims[3:4])*diff(flims[3:4]) + flims[3]
        
        # compute the bounds
        bounds  <- inset.size*diff(flims)[c(1,1,3,3)] + rep(c(figx, figy), each=2)        
        
        # set the background to the background colour
        if (!is.null(inset.bg)){
            b.xy    <- c(inset.x,inset.y,inset.x,inset.y) + inset.size[c(1,3,2,4)]*diff(lims)[c(1,3,1,3)]
            rect(b.xy[1], b.xy[2], b.xy[3], b.xy[4], border=NA, col=inset.bg)
        }
                
        par(plt=bounds, new=T)

        fun(...)
        
    } else {
        stop("Desired position not within plot region")
    }
}

#' @rdname inset
#' @param oldpar plotting parameters to which region is reset to
#' @param mfcol logical, should support for multiple panels be added?
#' @export
reset_pos <- function(oldpar, mfcol=TRUE){
    par(oldpar)
    if (prod(oldpar$mfrow) > 1){
      mfg <- oldpar$mfg
      mcol <- oldpar$mfcol
      if (mfcol){
        newmfg <- c(mfg[1] %% mcol[1] + 1 , (mfg[2] + mfg[1] %/% mcol[1] - 1) %% mcol[2] + 1)
      } else {
        newmfg <- c((mfg[1] + mfg[2] %/% mcol[2] - 1) %% mcol[1] + 1, mfg[2] %% mcol[2] + 1)
      }      
      if (!all(newmfg == c(1,1))) par('mfg' = newmfg)
    }
}


