#' Function to approximate Bayesian model probability using importance sampling.
#' The output logBML is approximating the Bayesian factor integral, and the support for one model over another, given the same data (i.e. the Bayes factor), can be approximated by exp(logBML[1] - logBML[2])
#' @param myf a Bayesian proportional probability density function
#' @param chain a chain from an MCMC fit. Exclude burning before submitting.
#' @param ndraws number of samples to use to approximate the bayesian model probability
#' @export
doBML <- function(myf,chain,ndraws=1e4){
  thetas <- apply(chain,2,mean);
  vcv    <- cov(chain);
  drws   <- rmvnorm(ndraws,mean=thetas,sigma=vcv);
  mbd  <- (sapply(1:dim(drws)[1],function(ii){myf(drws[ii,])}))
  dmv  <- dmvnorm(drws,mean=thetas,sigma=vcv,log=T)
  mod_hat <- myf(thetas)
  norm_hat <- dmvnorm(thetas,mean=thetas,sigma=vcv,log=T)
  logBML <- log(1/ndraws) + (-norm_hat+mod_hat) +
    log(sum(exp(mbd-mod_hat)/exp(dmv-norm_hat)))

  return(list(model_prob=mbd,normal_prob=dmv,draws=drws,thetas=thetas,vcv=vcv,logBML=logBML))
}
