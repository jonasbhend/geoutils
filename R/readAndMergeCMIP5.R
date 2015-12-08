#' Read in and collate a collection of NetCDF files.
#'
#' This function reads in data in CMIP5-like format and merges the
#' corresponding 'historical' and 'rcp' runs.
#'
#' @param files vector of pathnames to NetCDF files (old directory structure)
#' @param complete_ens logical, should models and runs be the same across RCPs?
#' @param complete_ts logical, should time series be complete?
#' @param mulc multiplicative factor (see \code{\link{readNetCDF}})
#' @param addc additive factor
#' @param valid_range range for quality control
#' @param validate logical, should quality control be perfomed (on a per model basis)?
#' @param experiments what experiment should be retained?
#' @param verbose logical, should progres notes be printed to screen?
#' @param startyear,endyear what years should time series be clipped to?
#' @param seas what seasons should be read in (see \code{\link{readSeasNetCDF}})
#' @param min.ncol Minimum number of time steps in historical to be included
#' @param varname variable name (see \code{\link{readNetCDF}})
#' @param mask logical, index of grid boxes to retain (see \code{\link{readNetCDF}})
#' @keywords utilities
#' @export
readAndMergeCMIP5 <- function(files, complete_ens=FALSE, complete_ts=TRUE, mulc=1, addc=0, valid_range=c(-1e10,1e10), validate=TRUE, experiments=c('historical', 'rcp26', 'rcp45', 'rcp60', 'rcp85'), verbose=TRUE, startyear=1900, endyear=2099, seas=c('djf', 'mam', 'jja', 'son'), min.ncol=length(seas), varname=NULL, mask=NULL){
  
  ## quick and dirty fix to get other indices up and running
  gseas <- seas
  if (any(!is.na(as.numeric(seas)))){
    gseas[!is.na(as.numeric(seas))] <- 'annual'
  }
  
  ## dissect file names
  filenames <- sapply(strsplit(files, '/'), function(x) x[length(x)])
  infotab <- data.frame(t(sapply(strsplit(gsub('.nc', '', filenames), '_'), function(x) x)))
  colnames(infotab) <- c('varname', 'mip_table', 'model', 'experiment', 'ensemble_member', 'grid')
  insti <- grep('IPCC', strsplit(files[1], '/')[[1]]) + 4
  infotab$institute <- sapply(strsplit(files, '/'), function(x) x[insti])
  infotab$filei <- seq(1,nrow(infotab))
  
  ## remove unwanted experiments
  infotab <- infotab[infotab$experiment %in% experiments,]
  
  ## work out joint constraints on files if globalmean is available
  grids <- unique(infotab$grid)
  if (sum(grids != 'globalmean') > 1) stop('Cannot deal with multiple grids')
  if (length(unique(infotab$varname[infotab$grid != 'globalmean'])) > 1) stop('Multiple variables')
  if (length(grids) > 1 & 'globalmean' %in% grids){
    splittab <- split(infotab, infotab[['grid']]) 
    globtab <- splittab[['globalmean']]
    infotab <- splittab[[grids[grids != 'globalmean']]]
    
    ## only use rows that have both global and regional data
    cnames <- c('model', 'experiment', 'ensemble_member')
    itxt <- apply(infotab[,cnames], 1, paste, collapse='')
    gtxt <- apply(globtab[,cnames], 1, paste, collapse='')
    infotab <- infotab[itxt %in% gtxt,]
    globtab <- globtab[gtxt %in% itxt,]    
  }
  
  ## generate variables for ease of use
  scenarios <- sort(unique(infotab$experiment[infotab$experiment %in% experiments]))
  if (is.null(varname)) varname <- unique(infotab$varname)
  models <- unique(infotab$model)
  
  ## apply constraints on joining models only for historical and rcp runs
  if ('historical' %in% scenarios){
    mtxt <- apply(infotab[,c('model', 'ensemble_member')], 1, paste, collapse='')
    unimtxt <- mtxt[mtxt %in% mtxt[infotab$experiment == 'historical'] & mtxt %in% mtxt[infotab$experiment != 'historical']]
    infotab <- infotab[mtxt %in% unimtxt,]
    if (length(grids) > 1 & 'globalmean' %in% grids){
      gtxt <- apply(globtab[,c('model', 'ensemble_member')], 1, paste, collapse='')
      globtab <- globtab[gtxt %in% unimtxt,]
    }
  }
  
  ## apply further constraints if all experiments have to be present for all models
  if (complete_ens){
    
    ## compare model/run combinations across scenarios
    mtxt <- apply(infotab[,c('model', 'ensemble_member')], 1, paste, collapse='')
    unimtxt <- unique(mtxt)
    for (scen in scenarios) unimtxt <- unique(mtxt[mtxt %in% unimtxt & mtxt %in% mtxt[infotab$experiment == scen]])
    if (length(unimtxt) == 0) stop('No complete sets of model data to read in')
    
    ## restrict tables
    infotab <- infotab[mtxt %in% unimtxt,]
    if (length(grids) > 1 & 'globalmean' %in% grids) {
      gtxt <- apply(globtab[,c('model', 'ensmeble_member')], 1, paste, collapse='')
      globtab <- globtab[gtxt %in% unimtxt,]
    }
  }
  
  ## initialise lists
  data <- list() # list to contain the regional area-averages
  if (length(grids) > 1 & 'globalmean' %in% grids) globdata <- list()
  ## loop over scenarios (not historical, as this is used for joins)
  for (scen in sort(scenarios[scenarios != 'historical'])){
    for (model in unique(infotab$model[infotab$experiment == scen])){
      for (runi in unique(infotab[['ensemble_member']][infotab$experiment == scen & infotab$model == model])){
        ## quick-n-dirty fix to not read in non-p1 runs
        if (length(grep('p1$', runi)) == 0) next
        if (verbose) print(paste(scen, model, runi))
        gtmp <- NULL ## initialise global temporary file for later checking
        ## read in scenario file for model, run, scenario
        dtmp <- try(readSeasNetCDF(files[infotab$filei[infotab$experiment == scen & infotab$model == model & infotab$ensemble_member == runi]], varname=varname, seas=tolower(seas), mulc=mulc, mask=mask) + addc)
        if (class(dtmp) == 'try-error') next
        if (ncol(dtmp) == min.ncol) next ## advance to next iteration immediately
        if (varname == 'hurs' & max(dtmp, na.rm=T) <= 1) dtmp <- dtmp*100
        if ("giorgi" %in% grids) dtmp <- select_region(dtmp, 1:40)
        if (substr(model, 1,4) == "GFDL" & varname == 'tas' & min(dtmp, na.rm=T) < 0) dtmp <- dtmp + 273.15
        if (diff(range(attr(dtmp, 'time'))) > 2000){
          print(paste('Time in', scen, model, runi, 'is wrong'))
          next ## advance to next iteration
        }
        if (length(grids) > 1 & 'globalmean' %in% grids) {
          gtmp <- try(readSeasNetCDF(files[globtab$filei[globtab$experiment == scen & globtab$model == model & globtab$ensemble_member == runi]], varname='tas', seas=tolower(gseas), mulc=1) + addc)
          if (class(gtmp) == 'try-error') next
          if (ncol(gtmp) == min.ncol) next ## advance to next iteration
          if (diff(range(attr(gtmp, 'time'))) > 2000){
            print(paste('Time in globalmean', scen, model, runi, 'is wrong'))
            next ## advance to next iteration
          }
        }
        ## if the scenario is split in two (e.g. historical and rcp85)
        if ('historical' %in% scenarios){
          hdtmp <- try(readSeasNetCDF(files[infotab$filei[infotab$experiment == 'historical' & infotab$model == model & infotab$ensemble_member == runi]], varname=varname, seas=tolower(seas), mulc=mulc, mask=mask) + addc)
          if (class(hdtmp) == 'try-error') next
          if (ncol(hdtmp) < min.ncol) next ## advance to next iteration immediately
          if (varname == 'hurs' & max(hdtmp, na.rm=T) <= 1) hdtmp <- hdtmp * 100
          if ("giorgi" %in% grids) hdtmp <- select_region(hdtmp, 1:40)
          if (substr(model, 1,4) == "GFDL" & varname == 'tas' & min(hdtmp, na.rm=T) < 0) hdtmp <- hdtmp + 273.15
          if (diff(range(attr(hdtmp, 'time'))) > 2000){
            print(paste('Time in historical', model, runi, 'is wrong'))
            next ## advance to next iteration
          }
          if (attr(hdtmp, 'units') != attr(dtmp, 'units')){
            print(paste('Units differ for', model, runi))
            next ## advance to next iteration
          }
          dtmp <- try(merge(hdtmp, dtmp), silent=TRUE)
          if (class(dtmp) == 'try-error') next
          rm(hdtmp)
          ## read in global mean temperature if needed
          if (length(grids) > 1 & 'globalmean' %in% grids) {
            hgtmp <- try(readSeasNetCDF(files[globtab$filei[globtab$experiment == 'historical' & globtab$model == model & globtab$ensemble_member == runi]], varname='tas', seas=tolower(gseas), mulc=1) + addc)
            if (class(hgtmp) == 'try-error') next
            if (ncol(hgtmp) == min.ncol) next ## advance to next iteration
            if (diff(range(attr(hgtmp, 'time'))) > 2000){
              print(paste('Time in globalmean historical', model, runi, 'is wrong'))
              next ## advance to next iteration
            }
            if (attr(gtmp, 'units') != attr(hgtmp, 'units')){
              print(paste('Units differ for', model, runi))
              next ## advance to next iteration
            }
            gtmp <- try(merge(hgtmp, gtmp), silent=TRUE)
            if (class(gtmp) == 'try-error') next
          }
        }
        
        if (complete_ts){
          inrange <- min(dtmp, na.rm=T) >= valid_range[1] & max(dtmp, na.rm=T) <= valid_range[2]
          xmask <- !apply(is.na(dtmp), 1, all)
          hasTime <- min(attr(dtmp, 'time')) <= startyear + 1 & max(attr(dtmp, 'time')) > endyear - 1
          hasMiss <- sum(is.na(dtmp[xmask,])) <= 5*length(seas)*sum(xmask)
          if (hasTime & hasMiss & inrange){
            data[[scen]][[model]][[runi]] <- dtmp
            if (length(grids) > 1 & 'globalmean' %in% grids) globdata[[scen]][[model]][[runi]] <- gtmp
          } else {
            if (verbose) {
              print(paste('Not read in:', model, runi, scen))
              if (!hasTime) print(paste('Time range:', paste(range(attr(dtmp, 'time'), na.rm=T), collapse='-')))
              if (!hasMiss) print(paste(sum(is.na(dtmp[xmask,]))/sum(xmask)/length(seas), 'years of missing values from', startyear,'to', endyear))
              if (!inrange) print(paste( model, runi, 'not in range:',paste(range(dtmp, na.rm=T), collapse=' to ')))
            }
          }
        } else {
          data[[scen]][[model]][[runi]] <- dtmp
          if (length(grids) > 1 & 'globalmean' %in% grids) globdata[[scen]][[model]][[runi]] <- gtmp
        }
        if (class(data[[scen]][[model]][[runi]]) != 'NetCDF') data[[scen]][[model]][[runi]] <- NULL
      }
    }        
  }
  
  ## diagnostic ouput
  if (verbose){
    print('finished reading in data')
    print(paste('Potential number of models:', length(unique(infotab$model))))
    for (i in seq(data)) print(paste(length(data[[i]]), 'models,', sum(sapply(data[[i]], length)),'runs in', names(data)[i]))
  }
  
  out <- list(data=data)
  if (length(grids) > 1 & 'globalmean' %in% grids) out$globdata <- globdata
  
  invisible(out)
  
}
