#' Function to extract effective sample sizes from a MCMC-CMR fit.
#' Calculated by using effectiveSize from the package coda.
#' @param fit A fit from MCMC_CMR. Extracts samples from fit$Chain
#' @return ESSs - an array of Effective sample size for each variable.
#' @export
ESS <- function(f1){
return(ESSs=coda::effectiveSize(coda::as.mcmc(f1$Chain[-c(1:dim(f1$Chain)[2]/2),])))
}
