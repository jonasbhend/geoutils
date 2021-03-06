#' @name select_models
#' @aliases select_ensembles select_ensmember select_ensmember2
#' 
#' @title
#' Select models and ensemble members
#' 
#' @description 
#' Functions to select models \code{select_models} and ensemble members
#' \code{select_ensembles} from an array of model results. With the corresponding
#' functions \code{select_ensmember} and \code{select_ensmember2} for backwards
#' compatibility.
#' 
#' @param x object of class 'NetCDF'
#' @param modi indices of models to select (first dimension in x, see details)
#' 
#' @details
#' The model selector can be either a vector of names corresponding to the row
#' names in \code{x}, a logical vector of length \code{nrow(x)}, or a numerical
#' vector with indices to retain.
#' 
#' @keywords utilities
#' @export
select_models <- function(x, modi=1:nrow(x)){
  ## convert the model denominator
  if (is.character(modi)){
    if (!is.null(rownames(x))){
      modi <- match(modi, rownames(x))
    } else {
      stop('Cannot select models by name - object has not row names')
    }
  } else if (is.logical(modi)){
    if (length(modi) != nrow(x)){
      stop('Cannot select models - logical vector of indices not of correct length')
    } else {
      modi <- which(modi)
    }
  } else if (is.numeric(modi)){
    if (max(modi) > nrow(x)) stop("Can't select the ensemble members you want")    
  }
  xtmp <- collapse2mat(x, first=TRUE)
  out <- array(xtmp[modi,], c(length(modi), dim(x)[-1]))
  attnames <- setdiff(names(attributes(x)), c('runs','dim', 'dimnames'))
  for (atn in attnames){
    attr(out, atn) <- attr(x, atn)
  }
  if (!is.null(attr(x, 'dimnames'))){
    adimnames <- attr(x, 'dimnames')
    if (!is.null(adimnames[[1]])) adimnames[[1]] <- adimnames[[1]][modi]
    attr(out, 'dimnames') <- adimnames
  }
  ## change run attribute if ensemble members are selected
  if (length(dim(x)) < 4 & !is.null(attr(x, 'runs'))){
    attr(out, 'runs') <- attr(x, 'runs')[modi]
  }
  invisible(out)
}
#' @rdname select_models
#' @param ensi indices (or text) of runs to select (second dimension of x)
#' @export
select_ensembles <- function(x, ensi=1:ncol(x)){
  xdims <- dim(x)
  xdims[2] <- length(ensi)
  if (is.numeric(ensi)){
    ## select ensemble members
    xtmp <- apply(x, seq(along=xdims)[-2], function(x) x[ensi])
    ## rearrange and write to array
    if (length(ensi) > 1){
      xout <- array(aperm(xtmp, c(2,1,3:length(dim(xtmp)))), xdims)
    } else {
      xout <- array(xtmp, xdims)
    }
  } else if (is.character(ensi)) {
    ## explicitly grab corresponding runs using the runnames attribute
    rnames <- attr(x, 'runnames')
    if (is.null(rnames)){
      stop('runnames not available')
    } else {
      xout <- array(NA, c(xdims[1:2], prod(xdims[-(1:2)])))
      xtmp <- array(x, c(dim(x)[1:2], prod(dim(x)[-(1:2)])))
      for (nn in 1:nrow(x)){
        for (ei in ensi){
          if (any(ei == rnames[nn,])){
            xout[nn,which(ensi == ei),] <- xtmp[nn,ei==rnames[nn,],]
          }
        }
      }
      xout <- array(xout, xdims)
      colnames(xout) <- ensi
    }
  } ## end of if on character vs numeric subscripts
  
  ## add attributes back in
  rownames(xout) <- rownames(x)
  atns <- attributes(x)[! names(attributes(x)) %in% c('dim', 'dimnames', 'runnames')]
  for (atn in names(atns)){
    attr(xout, atn) <- atns[[atn]]
  }
  invisible(xout)
}
#' @rdname select_models
#' @export
select_ensmember <- select_models
#' @rdname select_models
#' @export
select_ensmember2 <- select_ensembles

