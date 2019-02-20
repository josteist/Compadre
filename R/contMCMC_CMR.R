#' Function to continue sampling a fit, without adaptation.
#'
#' After an initial trial with \link{MCMC_CMR} which has seemingly reasonable mixing, this function
#' can be used to continue sampling using the estimated covariance proposal matrix.
#'
#' @param fit is the output from \link{MCMC_CMR}
#' @param x0 is the initial parameter values for continued sampling. If not given taken as the last parameter set from \emph{fit}
#' @param niter is the number of iterations
#' @param nthin is the thinning of the chain
#' @return a fit structure of class \emph{CMR_fit}
#' @export

# Making a 'continue sampling' function for MCMCs

contMCMC_CMR <- function(fit,x0=NULL,niter=1e5,nthin=10,...){
  # Just use MCMC_CMR with adapt=F, extract x vals and output
  # only new samples.
  # Can input other x0, for multi-chain runs, e.g.
  if (is.null(x0)){
  x0 <- fit$Chain[dim(fit$Chain)[1],];
  }

  cvstp <- fit$Covs;

  # How to extract settings from a 'call' as string;
  # We know that the call has MCMC_CMR( args...)
  # so we can remove the first 9 chars and the last, then split by ,
  # I don't think we need any of them. vmin is incorporatedin the covs and will not be adapted anyways.
  # tmpstr <- strsplit(substr(fit1$Call,start=10,stop=nchar(fit1$Call)-1),",")

  # which(grepl('niter',tmpstr[[1]]))

  # as.numeric(strsplit(tmpstr[[1]][2],'=')[[1]][2])

  newout <-   MCMC_CMR(cmrModel=fit$Model,niter=niter,nthin=nthin,
                       x0=x0,
                       cvstp = cvstp,
                       adapt=FALSE,...)
  return(newout)
}
