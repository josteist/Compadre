#' @export
print.CMR_fit <- function(fit){
  cat('A model fit using Compadre.\n')
  cat('\n')
  print.CMR_model(fit$Model)
  cat('\n')
  cat('  == Call to MCMC_CMR  == \n')
  print(fit$Call)
  cat('\n')
  cat('Size of chain:',dim(fit$Chain),'\n')
  cat('Type plot(fit) to plot estimated rates.')
  cat('Type hist(ESS(fit)) to plot histogram of effective sample sizes.')

}
# Making a 'name' for the model
# A CMR model with length(dts) intervals and dim(Obs)[1] observed taxa.
# Model has
