#' Function to extract effective sample sizes from a MCMC-CMR fit.
#' Calculated by using effectiveSize from the package coda.
#' @param fit A fit from MCMC_CMR. Extracts samples from fit$Chain
#' @return ESSs - an array of Effective sample size for each variable.
#' @export
ESS <- function(cmrfit){
return(ESSs=coda::effectiveSize(coda::as.mcmc(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1])/2),])))
}
