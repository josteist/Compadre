#' Function to approximate Bayesian model probability using importance sampling.
#' The output logBML is approximating the Bayesian factor integral, and the support for one model over another, given the same data (i.e. the Bayes factor), can be approximated by exp(logBML[1] - logBML[2])
#' @param cmrfit A fit of a model using Compadre. The function extracts the Bayesian probability density function and uses the last half of the chain.
#' @param ndraws number of samples to use to approximate the Bayesian model probability
#' @export
doBML <- function(cmrfit,ndraws=1e4){
  myf <- cmrfit$Model$probfun
  chain <- cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),]
  thetas <- apply(chain,2,mean);
  vcv    <- cov(chain);
  # mvtnorm::rmvnorm(draweps,sigma=cvstp)
  drws   <-  mvtnorm::rmvnorm(ndraws,mean=thetas,sigma=vcv);
  mbd  <- (sapply(1:dim(drws)[1],function(ii){myf(drws[ii,])}))
  dmv  <-  mvtnorm::dmvnorm(drws,mean=thetas,sigma=vcv,log=T)
  mod_hat <- myf(thetas)
  norm_hat <-  mvtnorm::dmvnorm(thetas,mean=thetas,sigma=vcv,log=T)
  logBML <- log(1/ndraws) + (-norm_hat+mod_hat) +
    log(sum(exp(mbd-mod_hat)/exp(dmv-norm_hat)))

  return(list(model_prob=mbd,normal_prob=dmv,draws=drws,thetas=thetas,vcv=vcv,logBML=logBML))
}
