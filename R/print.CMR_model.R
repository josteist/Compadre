#' @export
print.CMR_model <- function(mod){
  if (is.null(mod$clade1inx)){
    cat('A Capture-Mark-Recapture model using ',dim(mod$Obs)[1],' taxa,')
    cat(' spanning ',dim(mod$Obs)[2],' intervals. \n')
    cat('The model was generated',mod$date,'\n')
    cat('and has', mod$npar,'parameters.\n')
    cat('\n')

    # If only 1 clade model
    cat(' ==  Model has',length(unlist(mod$aix)),'drivers. == \n')
    if (length(unlist(mod$aix))>0){
      cat('Speciation rate: ',length(mod$aix$Covar_Speciation),'external driver(s)\n')
      cat('Extinction rate: ',length(mod$aix$Covar_Extinction),'external driver(s)\n')
      cat('Sampling rate:   ',length(mod$aix$Covar_Sampling),'external driver(s)\n')
      cat('Diversity dependence is',ifelse(length(mod$aix$DivDep_Speciation)>0,'ON','OFF'),'for speciation ')
      cat('and',ifelse(length(mod$aix$DivDep_Extinction)>0,'ON','OFF'),'for extinction rates.')
    }
    # cat('\n')
    if (sum(sapply(mod$reix,length)>0)==3){
      cat('All rates vary over time.')
    } else if (sum(sapply(mod$reix,length)>0)==2){
      cat(c('Speciation','Extinction','Sampling')[which(sapply(mod$reix,length)>0)[1]],' and ',
          c('speciation','extinction','sampling')[which(sapply(mod$reix,length)>0)[2]],' rates vary over time.')
    } else if (sum(sapply(mod$reix,length)>0)==1){
      cat(c('Speciation','Extinction','Sampling')[which.min(sapply(mod$reix,length)>0)],' rates vary over time.')
    } else {
      cat('\n')
      cat('No rates vary over time (i.e. have random terms).')
    }
  } else {
    cat(c('This is a model of interacting clades'))
  }

}

# Making a 'name' for the model
# A CMR model with length(dts) intervals and dim(Obs)[1] observed taxa.
# Model has
