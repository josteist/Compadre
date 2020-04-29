#' An adaptive Markov chain Monte Carlo sampler for estimating parameters of a CMR_model.
#'
#' After generating a model using \link{make_BayesCMR}, this function will sample the posterior. Default settings work relatively well for smaller models, but output needs to be checked for convergence (and more iterations are often necessary for models with many parameters). The sampler uses a covariance scaling to achieve good mixing. First half of the iterations are burning, and continually the proposals widths are tuned. The covariance structure of the samples in this burning phase is used for proposals in the latter half, insert ref here.
#'
#' @param cmrModel is a model generated with make_BayesCMR
#' @param niter is number of iterations. Defaults to 1e3.
#' @param nthin sets the thinning, i.e. samples for each nthin iterations are stored as output. Defaults to 10.
#' @param vmin sets the minimum and initial standard deviation of the normal proposals.
#' @param draweps number of iterations between each update of the proposal covariance matrix. Defaults to niter/100
#' @param x_init Potential initial values for the chain. If not provided the initial mean rates are optimes using the probability density function and all other parameters set to 0.
#' @param cvstp Covariance structure for proposal. This is tuned during burnin. Defaults to (2.38/(sqrt(cmrModel$npar)))^2*diag(vmin,cmrModel$npar))
#' @param adapt Adapt the covariance proposals during first half?
#' @return a fit structure with $Chain for samples, $Probs for posterior probabilities, $Accept number of accepted proposals in each block, $Model is the CMRmodel supplied as input, $Covs is the proposal covariance structure used in the last half of the chain.
#' @export
#'
MCMC_CMR <- function(cmrModel,niter=1e4,nthin=10,vmin=1e-5,
                         draweps=niter/noblc,noblc = 100,
                         x0=NULL,
                         cvstp = (2.38/(sqrt(cmrModel$npar)))^2*diag(vmin,cmrModel$npar),
                         adapt=TRUE){
  # Perhaps this works better if we ignore the 'covariance' part in the early phase (since most is prob directional anyways) and introduce it
  # Adaptive MCMC approach, tuning the stps to the covariance of the
  # chain for the first half of the run.
  # The covariance proposal is updated each draweps iteration the
  # first half og the total niter [number of iterations.]
  # _v2 should store all samples for each block(draweps), and use the samples from the last block
  # only to update the cov. Then we need a temp_X which stores ALL samples for each block.
  if (is.null(cmrModel$Clade1Mod)){
    # If this is null, then the model is for one clade only.
    if (is.null(x0)){
      # If no initial given, then optimize for main three parameters and use those as inits + only 0's
      # c(-2,-2.2,-2.1,rep(0,cmrModel$npar-3))
      x0 = c(optim(c(-1,-1.1,-1.2),function(x){-cmrModel$probfun(c(x[1],x[2],x[3],rep(0,cmrModel$npar-3)))})$par,
             rep(0,cmrModel$npar-3))
    }
  } else {
    # This is a two clade model.
    if (is.null(x0)){
      # If no initial given, then optimize for main three parameters and use those as inits + only 0's
      # c(-2,-2.2,-2.1,rep(0,cmrModel$npar-3))
      # Since we now can have drivers interacting the individual models wont work.
      tmpm1 <- make_BayesCMR(cmrModel$Clade1Mod$Obs,cmrModel$Clade1Mod$dts)
      tmpm2 <- make_BayesCMR(cmrModel$Clade2Mod$Obs,cmrModel$Clade2Mod$dts)
      xtmp1<-optim(c(-1.1,-1.2,-1),function(x){-tmpm1$probfun(x)})$par
      xtmp2<-optim(c(-1.1,-1.2,-1),function(x){-tmpm2$probfun(x)})$par #-cmrModel$Clade2Mod$probfun(c(x[1],x[2],x[3],rep(0,cmrModel$Clade2Mod$npar-3)))})$par

      x0 = c(xtmp1,rep(0,cmrModel$Clade1Mod$npar-3),
             xtmp2,rep(0,cmrModel$Clade2Mod$npar-3),
             rep(0,cmrModel$npar-cmrModel$Clade1Mod$npar - cmrModel$Clade2Mod$npar))
    }
  }
  fullcall <- deparse(match.call());
  # print(fullcall)
  print(paste0('MCMC run started at ', Sys.time()))
  x=x0;
  # noblc = niter/draweps;
  # Number of blocks; the whole chain is divided into blocks for faster computation. The proposal
  # var/covariance matrix is updated in the first half of the MCMC chain at the end of each block.
  stp <- niter/nthin # NOT USED?

  p_now = cmrModel$probfun(x) # Model probability now.
  X_Big = array(NA,c(ceiling(niter/nthin),cmrModel$npar)) # Large array for storing the parameter chains
  P_Big = array(NA,c(ceiling(niter/nthin))); # Long array to store the model probabilities.
  A_Big = array(NA,c(ceiling(niter/nthin))); # Long array to store the acceptance rates
  X_Big[1,]=x; # Inputting the initial parameters
  P_Big[1] = p_now; # and the initial model probability
  tmp_X = array(NA,c(draweps,cmrModel$npar)); # Temporary array to store all samples here for each block.
  # This is used to update the var/covar proposal matrix.

  Acc = 0; # Counter for acceptance of proposals
  prgb<- txtProgressBar(min=0,max=1,initial=0,style=3) # Print a progress bar.
  ii = 1; # counter for all iterations. Used to update the progres bar.
  tix = 1; # counter for the stored samples, i.e. for each nthin sample
  # adapt = T; # TRUE/FALSE if the var/covar should be updated.
  tmpA = 0;  # Temporary array to cound acceptances.
  # PErhasp rather than a vmin ,have a v_init = 100*vmin? OR the other way around, just initiaize with vmin*100`?`
  for (jj in 1:noblc){

    if (jj>1 & adapt==T){
      # Should augment this to be only
      # cvstp <- (2.38/(sqrt(cmrModel$npar)))^2*(cov(X_Big[(tix-draweps/nthin):(tix-1),])  +
      # diag(vmin,cmrModel$npar))
      # 261018. Below the new cov is based on the last 'half' of the stored
      # chain. This could be 'less', i.e. only using the cov based on the last
      # fewer number of iterations. I guess it should actually be based on all
      # samples the last n iterations, and not only the stored ones?
      # old:
      # cvstp <- (2.38/(sqrt(cmrModel$npar)))^2*(cov(X_Big[round(tix/2):(tix-1),])  +
      # diag(vmin,cmrModel$npar))
      # BELOW WAS COMMENTED OUT 110619 for testing setting the stepsize for main pars to /10 of others
      # cvstp <- (2.38/(sqrt(cmrModel$npar)))^2*(cov(tmp_X)  +
                                                 # diag(vmin,cmrModel$npar))
      # Testing 301019, instead of +vmin, do pmax(vmin)
      cvstp <- (2.38/(sqrt(cmrModel$npar)))^2*(cov(tmp_X)) + diag(vmin,cmrModel$npar);
      # diag(cvstp) <- pmax(diag(cvstp),vmin)


      # print(cvstp[1:4,1:4])
      if (jj>(noblc/2)){
        adapt=F;
      }

    }
    epss <- mvtnorm::rmvnorm(draweps,sigma=cvstp)
    for (ll in 1:draweps){
      ii = ii+1;
      p_new = cmrModel$probfun(x+epss[ll,])
      if (runif(1)<exp(p_new-p_now)){
        x = x+epss[ll,];
        p_now = p_new;
        Acc = Acc+1;
        tmpA = tmpA+1;
      }
      tmp_X[ll,] = x;
      if ((ii %% nthin)==0){
        X_Big[tix,]=x;
        P_Big[tix] = p_now;
        A_Big[tix] = tmpA;
        tmpA = 0;
        tix =tix+1;
        setTxtProgressBar(prgb,ii/niter);
      }
      # show(p_now)
    }

  }
  print(paste0('MCMC run finished at ', Sys.time()))
  out <- list(Call = fullcall,Chain=X_Big,Probs=P_Big,Accept=A_Big,Model=cmrModel,Covs = cvstp,date=date())
  attr(out,"class") <- "CMR_fit";
  return(out)
}
