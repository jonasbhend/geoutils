#' @name condense_ipcc
#' @aliases condense_ipcc
#' 
#' @title condense a list of NetCDF objects into an array
#' 
#' @description
#' These functions condense and subset a list of NetCDF objects (each containing time 
#' series of one model, ensemble run) into an array with the first two dimensions
#' being the model and runs respectively.
#' 
#' @param x list with objects of type NetCDF
#' @param year.range range of years for output data
#' @param seed random seed if missing values at end should be infilled
#' @param expand.ts logical, should missing values at end be infilled?
#' 
#' @keywords utilities
#' @export
condense_ipcc <- function(x, year.range=c(1870, 2100), seed=11, expand.ts=FALSE) {
  
  if (is.list(x) & all(sapply(x, function(y) 'time' %in% names(attributes(y))))){
    
    years <- sort(unique(unlist(lapply(x, attr, 'time'))))
    nseas <- median(table(floor(years)))
    out.arr <- array(NA, c(length(x), apply(sapply(x, dim), 1, max)[-nrow(sapply(x, dim))], length(years)))
    if (!is.null(names(x))) {
      xnames <- names(x) 
      rownames(out.arr) <- xnames[sort(tolower(xnames), ind=TRUE)$ix]
      if (any(xnames == 'r1i1p1')){
        runnum <- as.numeric(gsub('i.*', '', gsub('r', '', xnames)))
        physnum <- as.numeric(gsub('.*p', '', xnames))
        rownames(out.arr) <- xnames[sort(runnum + physnum/max(physnum), ind=T)$ix]
        xnames <- rownames(out.arr)
      }
    } else {
      xnames <- as.charcter(seq(along=x))
      rownames(out.arr) <- xnames
      names(x) <- xnames
    }
    ## write to common array
    x.dims <- sapply(x, dim)
    year.list <- lapply(x, attr, 'time')
    
    if (length(dim(out.arr)) == 4){
      
      runnames <- if (all(x.dims[1,] == 1)) array(NA, c(ncol(x.dims), max(sapply(x, function(y) length(attr(y, 'runnames')))))) else as.matrix(out.arr[,,1,1])
      rownames(runnames) <- rownames(out.arr)
      for (nn in xnames){
        out.arr[nn,1:x.dims[1,nn], 1:x.dims[2,nn], years %in% year.list[[nn]]] <- x[[nn]]
        if (all(x.dims[1,] == 1)){
          runnames[nn, 1:length(attr(x[[nn]], 'runnames'))] <- attr(x[[nn]], 'runnames')
        } else {
          runnames[nn,1:x.dims[1,nn]] <- attr(x[[nn]], 'runnames')
        }
      }
      out.arr <- out.arr[,,,years >= year.range[1] & years < year.range[2] + 1,drop=FALSE]
      
      ## constrain to year range (1870 - 2100)
      years <- years[years >= year.range[1] & years < year.range[2] + 1]
      ## mask out simulations with missing values to condense array
      mna <- apply(!is.na(out.arr[,,,years > 1900, drop=F]), 1:3, mean)
      mna[mna == 0] <- NA
      mthresh <- 1 - 1.5*(1 - median(mna, na.rm=T))
      ## mask out simulations with less than mthresh non-missing values
      for (i in seq(along=x)){
        xmask <- !apply(is.na(out.arr[i,,,,drop=FALSE]), 3, all)
        isim <- which(apply(!is.na(out.arr[i,,xmask,years > 1900,drop=FALSE]), 2, sum)/length(out.arr[i,1,xmask,years > 1900]) > mthresh)
        if (length(isim) != 0){
          out.arr[i,setdiff(1:ncol(out.arr),isim),,] <- NA
          out.arr[i,1:length(isim),,] <- out.arr[i,isim,,]
          if (!all(x.dims[1,] == 1)){
            runnames[i,setdiff(1:ncol(out.arr), isim)] <- NA
            runnames[i,1:length(isim)] <- runnames[i,isim]
          }
        }
      }
      ## reduce data array
      out.arr <- out.arr[,apply(!is.na(out.arr), 2, any),,,drop=FALSE]
      runnames[is.na(runnames)] <- ''
      attr(out.arr, 'runnames') <- runnames
    } else {
      for (nn in xnames){
        out.arr[nn,1:x.dims[1,nn], years %in% year.list[[nn]]] <- x[[nn]]
      }
      out.arr <- out.arr[,,years >= year.range[1] & years < year.range[2] + 1,drop=FALSE]
      ## constrain to year range (1870 - 2100)
      years <- years[years >= year.range[1] & years < year.range[2] + 1]
      attr(out.arr, 'runnames') <- rownames(out.arr)
    }
    
    
    if (expand.ts & length(dim(out.arr)) == 4){
      # extend every simulation to 2100 using (randomized) linear extrapolation
      # set seed in order to get comparable results
      set.seed(seed)
      for (i in seq(along=x)){
        nruns <- sum(apply(!is.na(out.arr[i,,,]), 1, any))
        xmask <- !apply(is.na(out.arr[i,,,]), 2, all)
        for (j in 1:nruns){
          datmax <- max(which(!is.na(out.arr[i,j,min(which(xmask)),])))
          if (datmax < dim(out.arr)[4]){
            usedims <- max(1, dim(out.arr)[4] - nseas*30+1)
            tmp <- out.arr[i,j,,usedims:dim(out.arr)[4]]
            tmp <- array(rep(tmp, each=2), c(2, nrow(tmp)*nseas, ncol(tmp)/nseas))
            tmp[2,,] <- rep(years[usedims:dim(out.arr)[4]], each=dim(out.arr)[3])
            tmpout <- apply(tmp, 2, function(x) {
              x.lm <- lm(x[1,] ~ x[2,])
              x.sd <- sd(x.lm$resid)
              x.pred <- coef(x.lm) %*% rbind(1, x[2,])
              xout <- x[1,]
              xout[is.na(xout)] <- rnorm(sum(is.na(xout)), x.pred[is.na(xout)], x.sd)
              return(xout)
            })
            out.arr[i,j,,usedims:dim(out.arr)[4]] <- t(tmpout)
          }
        }
      }
    }
    attnames <- setdiff(names(attributes(x[[1]])), c('dim', 'dimnames', 'runnames'))
    for (atn in attnames){
      attr(out.arr, atn) <- attr(x[[1]], atn)
    }
    if (!is.null(attr(x[[1]], 'nens'))) attr(out.arr, 'nens') <- sapply(x, attr, 'nens')
    if (!is.null(attr(x[[1]], 'runs'))) attr(out.arr, 'runs') <- lapply(x, attr, 'runs')
    attr(out.arr, 'time') <- years
    return(out.arr)
  } else {
    warning('condense_ipcc does not work with this object')
    return(x)
  }
}

#' @rdname condense_ipcc
#' @export
condense_picntrl <- function(x) {
  x.dims <- sapply(x, dim)
  out.arr <- array(NA, c(length(x), apply(x.dims, 1, max)))
  times <- unique(sort(unlist(lapply(x, attr, 'time'))))
  nseas <- 12/median(table(floor(times)))
  times <- times[(times %% 1) %% (nseas/12) == 0]
  outtimes <- seq(times[1], length=dim(out.arr)[4], by=nseas/12)
  for (i in seq(along=x)){
    # check whether first time step corresponds to first timestep in list
    time.offset <- min(which(attr(x[[i]], 'time') %% 1 == outtimes[1] %% 1))
    out.arr[i,1:x.dims[1,i], 1:x.dims[2,i], seq(time.offset,x.dims[3,i])-time.offset + 1] <- x[[i]][,,time.offset:x.dims[3,i]]
  }
  # reduce data array
  timmask <- apply(!is.na(out.arr), 4, any)
  out.arr <- out.arr[,apply(!is.na(out.arr), 2, any),,timmask]
  attr(out.arr, 'time') <- outtimes[timmask]
  return(out.arr)
}
#' @rdname condense_ipcc
#' @param x list of objects of type NetCDF
#' @param clean logical, should resulting array be cleaned up?
#' @export
implode <- function(x, nfrac=0.8, clean=FALSE) {
  if (!is.list(x)) stop('implode works on list objects only')
  times <- sort(unique(unlist(lapply(x, function(y) attr(y, 'time')))))
  xdims <- sapply(x, function(y) dim(y))
  maxdim <- apply(xdims, 1, max)
  out <- array(NA, c(length(x), maxdim[seq(1,length(maxdim)-1)], length(times)))
  if (length(maxdim) == 2){
    for (i in seq(along=x)){
      out[i,1:xdims[1,i],times %in% attr(x[[i]], 'time')] <- x[[i]]
    }
  } else if (length(maxdim) == 3){
    for (i in seq(along=x)){
      out[i,1:xdims[1,i],1:xdims[2,i],times %in% attr(x[[i]], 'time')] <- x[[i]]
    }
  }
  rownames(out) <- names(x)
  ## remove runs with missing data
  if (clean){
    if (length(maxdim) == 2){
      out <- out[apply(!is.na(out), 1, mean) > nfrac,,,drop=F]
    } else if (length(maxdim) == 3) {
      out <- out[apply(!is.na(out), 1, mean) > nfrac,,,,drop=F]
    }
  }
  attnames <- setdiff(names(attributes(x[[1]])), c('dim', 'dimnames'))
  for (atn in attnames){
    attr(out, atn) <- attr(x[[1]], atn)
  }
  attr(out, 'time') <- times
  return(out)
}
