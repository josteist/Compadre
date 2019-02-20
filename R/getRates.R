#' Extract macroevolutionary rates from CMR model fit.
#'
#' @param fit A fit returned from \link{MCMC_CMR}
#' @param ns number of samples used in chain.
#' @export
getRates <- function(fit,ns=250){
  # repar this so as to also take in list of fits?
  tpl <- ceiling(seq(dim(fit$Chain)[1]/2,dim(fit$Chain)[1],length.out=ns));
  # extracting spc rate quantile
  SpecRates <-exp((sapply(1:length(tpl),function(ii){
    fit$Model$ratefunc[[1]](fit$Chain[tpl[ii],])+
      fit$Model$ratefunc[[4]](fit$Chain[tpl[ii],])})))
  # extinction quantile
  ExtRates <- exp((sapply(1:length(tpl),function(ii){
    fit$Model$ratefunc[[2]](fit$Chain[tpl[ii],])+
      fit$Model$ratefunc[[5]](fit$Chain[tpl[ii],])})))

  # sampling quantiles
  SampRates <- exp((sapply(1:length(tpl),function(ii){
    fit$Model$ratefunc[[3]](fit$Chain[tpl[ii],])})))
  return(list(SpecRates = SpecRates, ExtRates=ExtRates,SampRates=SampRates))
}
