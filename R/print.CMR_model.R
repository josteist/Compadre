#' @export
print.CMR_model <- function(mod){
  if (is.null(mod$Clade1Mod)){
    cat(' == Compadre model == \n')
    cat('Model includes ',dim(mod$Obs)[1],' taxa,')
    cat(' spanning ',dim(mod$Obs)[2],' intervals. \n')
    cat('The model was generated',mod$date)
    cat(' and has', mod$npar,'parameters.\n')
    cat('\n')

    # If only 1 clade model
    cat(' ==  Model terms  == \n')
      cat('Speciation',as.character(mod$spec),'\n')
      cat('Extinction',as.character(mod$ext),'\n')
      cat('Sampling',as.character(mod$samp),'\n')


    cat('\n')
  } else {
    cat(' == Compadre model of interacting clades == \n')
    cat('Model includes',dim(mod$Clade1Mod$Obs)[1],'taxa in clade 1 and',dim(mod$Clade2Mod$Obs)[1],'in clade 2,')
    cat(' spanning ',dim(mod$Clade1Mod$Obs)[2],' intervals. \n')
    cat('The model was generated',mod$date)
    cat(' and has', mod$npar,'parameters.\n')
    cat('\n')

    # If only 1 clade model
    cat(' ==  Model terms  == \n')
    cat('Speciation clade 1',as.character(mod$Clade1Mod$spec),'\n')
    cat('Extinction clade 1',as.character(mod$Clade1Mod$ext),'\n')
    cat('Sampling clade 1',as.character(mod$Clade1Mod$samp),'\n')
    cat('\n')
    cat('Speciation clade 2',as.character(mod$Clade2Mod$spec),'\n')
    cat('Extinction clade 2',as.character(mod$Clade2Mod$ext),'\n')
    cat('Sampling clade 2',as.character(mod$Clade2Mod$samp),'\n')

      }

}

# Making a 'name' for the model
# A CMR model with length(dts) intervals and dim(Obs)[1] observed taxa.
# Model has
