
#' Function to construct CMR-model.
#'
#' This function generates a structure containing the necessary functions for a CMR analysis of a fossil dataset.
#'
#' @param Obs a matrix with size \emph{number of taxa} by \emph{number of intervals (n)}. Each taxa has a row with 0's (unobserved) and 1's (observed) for each interval in the analysis. Oldest interval is first column.
#' @param dts a vector of interval durations. Defaults to a series of 1's if not supplied. Oldest interval is first entry.
#' @param RE  a vector of TRUE/FALSE statements whether or not to allow rates to vary across intervals for speciation, extinction and sampling rates, respectively. Defaults to c(FALSE, FALSE, FALSE)
#' @param SpecTS a vector or matrix of drivers for speciation rate. One row for each driver and one column for all but the last interval (i.e. \emph{n-1}, since there is no rate transitioning out of the last interval). Defaults to NULL. Oldest interval is first column.
#' @param ExtTS a vector or matrix of drivers for extinction rates. Defaults to NULL. Oldest interval is first column. Each driver of length \emph{n-1}.
#' @param SmpTS a vector or matrix of drivers for sampling rates. One row for each driver and one column for each interval. Defaults to NULL. Oldest interval is first column. Each driver has \emph{n} entries.
#' @param DivDep a switch to set inclusion of diversity dependence in speciation and extinction rates, respectively. Defaults to c(FALSE, FALSE)
#' @param pfix a switch to select which approach is used to solve the identifiability problem in the model. If \emph{pfix = 1}, then the sampling rate in the first and the last intervals are assumed to be equal to the mean sampling rate for the whole period. If \emph{pfix = 2}, then the first two intervals, and the last two intervals have the same sampling rate. Defaults to pfix = 2.
#' @param priorsNorm_Cov sets the parameters for the normal prior for the impact of the drivers. Defaults to Norm(mu=0,sd=2).
#' @param priorsNorm_Mus sets the parameters for the normal prior for the (log) mean speciation, extinction and sampling rates. Given as a list of three entries with (mu,sd). Defaults to rep(list(c(-4,4)),3), i.e. all normal priors with mu=-4 and sd=3.
#' @param priorsUnif sets the parameters for the uniform prior for the log variances of the random efSpecTS. Same for all three variances, dunif(-3,3).
#' @param replRE_1 a switch for how to deal with very special drivers. Not fully implemented or tested properly.
#' @return The function returns an object of class CMR_model. This most important output is out$probfun which is a function for the posterior. This output can be fed directly into the sampler \link{MCMC_CMR}. It also returns \emph{Obs} and \emph{dts} as well as the full call (exlcluding default settings), \emph{call}
#' @export
#' @examples mod1 <- make.BayesCMR(Obs)
#' fit <- MCMC_CMR(mod1)
#' matplot(fit$Chain[,1:3],type="l")
#' This version implements the speciation and extinction rates using the length of the bins, not the
#' inter-bin midpoint difference.

make_BayesCMRv2 <- function(Obs,dts=rep(1,dim(Obs)[2]),
                             RE=c(FALSE,FALSE,FALSE),
                             SpecTS=NULL,ExtTS=NULL,SmpTS=NULL,
                             DivDep=c(FALSE,FALSE),pfix=2,priorsNorm_Cov=c(0,2),
                             priorsNorm_Mus = rep(list(c(-4,4)),3),
                             priorsUnif=c(-3,3),replRE_1 = F){
  # Model generating function for a Compadre analysis.
  # Minimum input is a matrix of observed/unobserved of dimensions
  # taxa by temporal interval. dts is vector of interval
  # durations (default sets to 1).

 # v2 is as main version of 26.03.2019, BUT with one exception. For seniorities
  # and extinction probabilities in the pradel likelihood the time between
  # midpoints is used and not the proper interval durations.
    # Input checks to code in
  # check dimensions of Obs vs dts, obs vs SpecTS,exts,SmpTS
  #
  ## INPUT CHECKS
  fullcall <- deparse(match.call());

  if (length(dts) != dim(Obs)[2]){
    stop("Number of intervals in *dts* does not match dimension of Obs")
  }
  if (!is.null(SpecTS)){
    # if drivers for speciation rates are given, check the dimensionality
    if (NCOL(SpecTS)!=(dim(Obs)[2]-1)){
      stop("Size of speciation rate drivers does not match dimension of Obs. SpecTS should be of NCOL==(dim(Obs)[2]-1)")
    }
  }

  if (!is.null(ExtTS)){
    # if drivers for speciation rates are given, check the dimensionality
    if (NCOL(ExtTS)!=(dim(Obs)[2]-1)){
      stop("Size of extinction rate drivers does not match dimension of Obs. ExtTS should be of NCOL==(dim(Obs)[2]-1)")
    }
  }

  if (!is.null(SmpTS)){
    # if drivers for speciation rates are given, check the dimensionality
    if (NCOL(SmpTS)!=(dim(Obs)[2])){
      stop("Size of sampling rate drivers does not match dimension of Obs. ExtTS should be of NCOL==(dim(Obs)[2])")
    }
  }
  # Obs - matrix of observations (absence/presence) of dimensions (species by intervals)
  # dts - duration of intervals (can be named?)
  # RE  - use Random efSpecTS for [speciation, extinction, sampling] rates respectively T/F(default)
  # xxxts - supplied covariate time-series of dimensions (#ts by #intervals for sampling covariates and #ts by #intervals-1 for spec/ext
  # DivDep- switch for diversity dependence in [speciation, extinction] T/F(default)
  # pfix - method for fixing sampling rates in first and last bin. pfix = 1, rates are same and = mean, pfix=2; first=second, and second-last=last
  # priorsNorm - prior parameters for all non-variance parameters (mean[rates], all covariate efSpecTS and divdep)
  # priorsUnif - prior parameters for variance terms [4:6] for random efSpecTS dunif([1],[2])
  # replRE_1   - if T and a timeseries with only 1 entry of 1 is supplied the Random effect for this entry is removed.
  nprs1 <- priorsNorm_Cov; # prior pars. for drivers
  if (length(priorsNorm_Mus)==3){
    # then it's a list of length three, most likely one set of pars for each
    nprs2 <- priorsNorm_Mus; # prior pars. for mean rates.
  } else if (length(priorsNorm_Mus)==2){
    # then just two, and should be repeated
    nprs2 <- rep(list(priorsNorm_Mus),3)
  }


  #All normal priors are µ=0 and sd=10
  priun<- priorsUnif;# Obs is 1/0 observation matrix. dts is (possibly) a vector of time-durations of the intervals
  # RE is a TRUE/FALSE array for inclusion of random-efSpecTS (in practice makes rates slightly independent across time)
  # xxxts are possible time series of covariates for fec(speciation rate), ext (extinction rate) or sampling (smp)
  # If xxtx are supplied, they will be included as temporal covariates.
  # they should have dimensions (n_cov by dim(Obs)[1]-1) for f & e and
  # dim(Obs)[1] for smp. If RE are true random efSpecTS will be added for fec,
  # ext and samp respectively. Should be able to handle both.
  # DivDep is a true/false statement linking speciation and/or extinction rates to richness, estimated as n_obs[i]/p[i], normalized to have mean 0 and sd=1 each calculation
  # Possible extensions:
  # -  allow for splitting Obs into groups with effect for each group.
  # - Have a switch for fixing the sampling rates at the extremes; now first and last
  #   don't have a RE, which makes the model identifiable. Could also be first and last have
  #   same RE as second and second-last?
  #       - pfix = 1; - first and last log(p_samp) = µ[p_samp]
  #       - pfix = 2; - first and last log(p_samp) = second and second last p_samp
  #   --> Done.
  # -
  # Improvement: - remove the RE if the TS has one interval which will have a special
  # rate, e.g. for testing if there is elevated ext. rate for a specifc transition. Now
  # havine a TS with one 1 will also have a RE for this interval.
  # 1. identify which par and which interval.
  # 2. Then it might be easier to make index arrays outside the if statements.
  # This will only be an issue under RE=T and supplied TS. DONE
  undv <- make_unvd(Obs);
  ts = list(SpecTS,ExtTS,SmpTS)
  # Calcing these here, in case they are accidentaly changed in the globalEnv.
  u =undv$u;
  v = undv$v;
  n = undv$n;
  d = undv$d;

  # There must be a better way to do this...

  prior <- list(); #making all prior terms part of this list.
  prix = 1; # ticker for prior parts.
  # can I define the prior cumulatively? just add no.
  # dim(array) yields NULL, which is annnoying.
  # checking the dimension of the covariates for each par
  origSpecTS = SpecTS;
  origSmpTS = SmpTS;
  origExtTS = ExtTS;
  if (!is.null(SpecTS)){
    # if not null, is it 1 or more cvts?
    if (is.null(dim(SpecTS))){
      nSpecTS = 1;
      # Normalizing: 290119 - made internal to function
      SpecTS = normfun(origSpecTS);
    } else {
      nSpecTS = dim(SpecTS)[1];
      SpecTS   = t(apply(origSpecTS,1,normfun))
    }
  } else {nSpecTS=0}
  if (!is.null(ExtTS)){
    # if not null, is it 1 or more cvts?
    if (is.null(dim(ExtTS))){
      nExtTS = 1;
      ExtTS = normfun(origExtTS);
    } else {
      nExtTS = dim(ExtTS)[1];
      ExtTS  = t(apply(origExtTS,1,normfun));
    }
  } else {nExtTS=0;}
  if (!is.null(SmpTS)){
    # if not null, is it 1 or more cvts?
    if (is.null(dim(SmpTS))){
      nSmpTS = 1;
      SmpTS = normfun(origSmpTS);
    } else {
      nSmpTS = dim(SmpTS)[1];
      SmpTS = t(apply(origSmpTS,1,normfun))
    }
  } else {nSmpTS=0}
  ts = list(SpecTS,ExtTS,SmpTS) # storing normalized drivers.
  ts_orig = list(origSpecTS,origExtTS,origSmpTS)


  nhps <- 3+sum(RE*1)+sum(DivDep*1) + nSmpTS+nExtTS+nSpecTS; # number of mean+hyperparameters
  nTS  <- c(nSpecTS,nExtTS,nSmpTS);
  # ALL x-arrays endd with the RE's, so all indexes below here that are not the first 3+sum(RE) sohuld have sum(DivDep) added to them.
  # All these have priors dnorm(0,5), the Sd[RE] are log(sd) really.
  # [µ_fec,µ_ext,µsamp] &
  # (possibly)[st_RE_fec,st_RE_ext,st_RE_samp] &
  # possibly [alphaTS1,alphaTS2 etc]
  #
  alphinx= list(); # These are effect of timeseries
  alphinx[[1]] <- seq(3+sum(RE*1)+sum(DivDep*1)+1,length.out=nSpecTS)
  alphinx[[2]] <- seq(3+sum(RE*1)+sum(DivDep*1)+1+length(alphinx[[1]]),length.out=nExtTS)
  alphinx[[3]] <- seq(3+sum(RE*1)+sum(DivDep*1)+1+length(alphinx[[1]])+length(alphinx[[2]]),length.out=nSmpTS)
  alphinx[[4]] <- seq(3+sum(RE*1)+1,length.out=DivDep[1]*1)
  alphinx[[5]] <- seq(3+sum(RE*1)+DivDep[1]*1+1,length.out=DivDep[2]*1)

  names(alphinx) <- c('Covar_Speciation','Covar_Extinction','Covar_Sampling',
                      'DivDep_Speciation','DivDep_Extinction')

  lns = c(length(dts)-1,length(dts)-1,length(dts)); # no pars to feed Pradel_unvd func
  tix = nhps+1; # ticker for no of params. These are for RE's only. TS efSpecTS have been included in nhps
  reix = list(); # list of RE indexes into x for speciation, extinction and sampling
  # If simplest model then nhps = 3
  if (nhps==3){
    myf <- function(x){
      # Simple 3 rate model. Estimated on log x. Fecundity
      # is assumed to be poisson-like, i.e. actual speciation/fecundity for an interval
      # of duration dt is fec*dt.
      lfec = rep(x[1],length(dts)-1);
      lext = rep(x[2],length(dts)-1);
      lsmp = rep(x[3],length(dts));

      ll <- pradel_unvd_gam(
        ext = sapply(1:(length(dts)-1),function(ii){
          (exp(lext[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2)))-1))/
            (exp(lfec[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))))-exp(lext[ii]))}),
        gam <- sapply(1:(length(dts)-1),function(ii){
          (1-(exp(lext[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2)))-1))/
             (exp(lfec[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))))-exp(lext[ii])))/
            exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))}),
        p   = sapply(1:(length(dts)),function(ii){rate2prob(exp(lsmp[ii]),dts[ii])}),
        u,n,v,d);

      # For different priors for the three global rates
      ll$LogL <- ll$LogL + sum(sapply(1:3,function(ii){dnorm(x[ii],nprs2[[ii]][1],nprs2[[ii]][2],log=T)}))
      if (is.infinite(ll$LogL)){ll$LogL = -Inf};
      if (is.na(ll$LogL)){ll$LogL = -Inf};
      if (is.nan(ll$LogL)){ll$LogL = -Inf};
      return(ll$LogL)

    }
    myll <- function(x){
      # Simple 3 rate model. Estimated on log x. Fecundity
      # is assumed to be poisson-like, i.e. actual speciation/fecundity for an interval
      # of duration dt is fec*dt.
      lfec = rep(x[1],length(dts)-1);
      lext = rep(x[2],length(dts)-1);
      lsmp = rep(x[3],length(dts));

      ll <- pradel_unvd_gam(
        ext = sapply(1:(length(dts)-1),function(ii){
          (exp(lext[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2)))-1))/
            (exp(lfec[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))))-exp(lext[ii]))}),
        gam <- sapply(1:(length(dts)-1),function(ii){
          (1-(exp(lext[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2)))-1))/
             (exp(lfec[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))))-exp(lext[ii])))/
            exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))}),
        p   = sapply(1:(length(dts)),function(ii){rate2prob(exp(lsmp[ii]),dts[ii])}),
        u,n,v,d);
      if (is.infinite(ll$LogL)){ll$LogL = -Inf};
      if (is.na(ll$LogL)){ll$LogL = -Inf};
      if (is.nan(ll$LogL)){ll$LogL = -Inf};
      return(ll)

    }

    ltmpfun <- list(); # list of functions for fec, ext and samp
    ltmpfun[[1]] <- function(x){rep(x[1],length(dts)-1)}
    ltmpfun[[2]] <- function(x){rep(x[2],length(dts)-1)};
    ltmpfun[[3]] <- function(x){rep(x[3],length(dts))};
    ltmpfun[[4]] <- function(x){0};
    ltmpfun[[5]] <- function(x){0}; # for completeness of output only.
    ptmpfun <- list();
    n_est <- function(x){ (n[1:length(dts)-1]/
                             (rate2prob(exp(ltmpfun[[3]](x)),dts)[1:(length(dts)-1)]))}
    n_norm <- function(x){( n_est(x)-mean(n_est(x)))/sd(n_est(x))};

  } else {
    # more complicated model
    ltmpfun <- list(); # list of functions for fec, ext and samp
    ptmpfun <- list(); # list of functions for PRIORs. For random efSpecTS only.
    # if no prior function is added here, place
    ptmpfun[[1]] <- function(x){0}
    ptmpfun[[2]] <- function(x){0}
    ptmpfun[[3]] <- function(x){0}
    # First three are always µ[fec],µ[ext],µ[samp]. If ~time the next three
    # are either effect (if TS is supplied) or std of Random EfSpecTS.
    for (jj in 1:3){
      # for each par [fec, ext, samp]
      # nTS
      # if ((RE[jj]==FALSE) & (is.null(ts[[jj]]))){
      if ((RE[jj]==FALSE) & (nTS[jj]==0)){
        # no temporal dependency AND no ts supplied. simple
        # ltmpfun[[jj]] <- force(jj) function(x){force(jj);rep(x[jj],lns[jj])}
        ltmpfun[[jj]] <- eval(substitute(function(x){rep(x[j1],lns[j1])},list(j1=jj)))
        # the below didn't work, might be similar problems in more complex stuff below
        # function(x){rep(x[jj],lns[jj])};
        # This might fix itself when it's run as function.
        # eval(substitute(function(x){rep(x[j1],lns[j1])},list(j1=jj)))
        # I don't understand how I can use this, now this returns ltmpfuns which
        # only uses jj=3, even though it's indexed as lttmpfun[jjj!=3]
        # This seems to work, but is quite the mouthful.

      } else { # if some temporal dependency in this par.
        # Changing this to n's, should perhaps be able to get an 'empty' array for listapply
        if (nTS[jj]==0){
          # if no external timeseries is given.
          if (jj==1){
            # ltmpfun[[jj]] <- eval(substitute(function(x){rep(x[j1],lns[j1])},list(j1=jj)))
            #
            ltmpfun[[1]] <-   eval(substitute(function(x){x[1] +
                c(  x[seq(stix,length.out=(length(dts)-1))])},list(stix=tix)))
            ptmpfun[[1]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-1))],0,exp(x[3+cumsum(RE)[1]]),log=T))}, list(stix=tix)))
            reix[[jj]] <- c(seq(tix,length.out=lns[jj]))
            tix = tix+lns[jj];

            # don't really know why, but the eval(substitute actually works.)
          } else if (jj==3){
            if (pfix==1){ #sampling in first and last interval is equal to µ[samp]
              ltmpfun[[3]] <- eval(substitute(function(x){x[3] +
                  c(0,x[seq(stix,length.out=(length(dts)-2))],0)},list(stix=tix)))
              ptmpfun[[3]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-2))],0,exp(x[3+cumsum(RE)[3]]),log=T))}, list(stix=tix)))
              reix[[jj]] <- c(seq(tix,length.out=lns[jj]-2))
              tix = tix+lns[jj]-2;
            }  else if (pfix==2){ # sampling first and last is eq to second and second last.
              ltmpfun[[3]] <- eval(substitute(function(x){x[3] +
                  c(x[stix],x[seq(stix,length.out=(length(dts)-2))],x[stix+length(dts)-3])},list(stix=tix)))
              ptmpfun[[3]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-2))],0,exp(x[3+cumsum(RE)[3]]),log=T))}, list(stix=tix)))
              reix[[jj]] <- c(seq(tix,length.out=lns[jj]-2))
              tix = tix+lns[jj]-2;

            }
          } else if (jj==2){
            ltmpfun[[2]] <- eval(substitute(function(x){x[2] +
                c(  x[seq(stix, length.out=(length(dts)-1))])},list(stix=tix)))
            ptmpfun[[2]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-1))],0,exp(x[3+cumsum(RE)[2]]),log=T))}, list(stix=tix)))
            reix[[jj]] <- c(seq(tix,length.out=lns[jj]))
            tix = tix+lns[jj];
          }
        } else {
          # if timeseries is given for this. Assume that the TS has same length
          if (RE[jj]==FALSE){
            if (nTS[jj]==1){
              # if only one covariate
              # tix is not right here, the TS efSpecTS are earlier in X, imagine some RE's have been included, then tix is not the value of the TS effect
              ltmpfun[[jj]] = eval(substitute(function(x){x[j1] + ts[[j1]]*x[alphinx[[j1]]]},list(j1=jj)))
              # if no RE's then the effect of this TS is parameter 4
              # tix = tix+1;
            } else {
              # more covariates for same parameter.
              ltmpfun[[jj]] <- eval(substitute(function(x){x[j1] +
                  rowSums(sapply(1:dim(ts[[j1]])[1],function(ii){x[alphinx[[j1]][ii]]*ts[[j1]][ii,]}))},list(j1=jj)))
              # tix = tix+dim(j)[1]
            }
          } else {
            # If also Random efSpecTS. Only here will the problem with
            # one RE AND one TS effect (if testing only one interval) be a problem.
            #
            if (nTS[jj]==1){
              # (is.null(dim(ts[[jj]]))){
              # if only one covariate
              if (jj==1){
                #replRE_1
                if (sum(ts[[1]]!=0)==1 & replRE_1){
                  # if only a time-series with 1 interval special. Then this interval should have NO random effect.
                  # This is not implemented if given > 1 TS.
                  ltmpfun[[1]] <-   eval(substitute(function(x){
                    x[1] + ts[[1]]*x[alphinx[[j1]]] +
                      c(x[seq(stix,length.out=(which(ts[[1]]>0)-1))],0,
                        x[seq(stix+which(ts[[1]]>0)-1,length.out=(length(dts)-which(ts[[1]]>0)-1))])},list(stix=tix,j1=jj)))
                  # c(  x[seq(stix,                 length.out=(length(dts)-1))])},list(stix=tix,j1=jj)))
                  ptmpfun[[1]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-2))],0,exp(x[3+cumsum(RE)[1]]),log=T))}, list(stix=tix)))
                  reix[[jj]] <- c(seq(tix,length.out=lns[jj]-1))
                  tix = tix+lns[jj]-1;
                } else {
                  ltmpfun[[1]] <-   eval(substitute(function(x){
                    x[1] + ts[[j1]]*x[alphinx[[j1]]] +
                      c(  x[seq(stix,                 length.out=(length(dts)-1))])},list(stix=tix,j1=jj)))
                  ptmpfun[[1]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-1))],0,exp(x[3+cumsum(RE)[1]]),log=T))}, list(stix=tix)))
                  reix[[jj]] <- c(seq(tix,length.out=lns[jj]))
                  tix = tix+lns[jj];
                }
              } else if (jj==2){
                if (sum(ts[[2]]!=0)==1 & replRE_1){
                  ltmpfun[[2]] <-   eval(substitute(function(x){
                    x[2] + ts[[2]]*x[alphinx[[2]]] +
                      c(x[seq(stix,length.out=(which(ts[[2]]>0)-1))],0,
                        x[seq(stix+which(ts[[2]]>0)-1,length.out=(length(dts)-which(ts[[2]]>0)-1))])},list(stix=tix,j1=jj)))
                  # c(  x[seq(stix,                 length.out=(length(dts)-1))])},list(stix=tix,j1=jj)))
                  ptmpfun[[2]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-2))],0,exp(x[3+cumsum(RE)[2]]),log=T))}, list(stix=tix)))
                  reix[[jj]] <- c(seq(tix,length.out=lns[jj]))
                  tix = tix+lns[jj]-1;

                } else {
                  ltmpfun[[2]] <- eval(substitute(function(x){x[2] + ts[[j1]]*x[alphinx[[j1]]]+ c(  x[seq(stix, length.out=(length(dts)-1))])},list(j1=jj,stix=tix)))
                  ptmpfun[[2]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-1))],0,exp(x[3+cumsum(RE)[2]]),log=T))}, list(stix=tix)))
                  reix[[jj]] <- c(seq(tix,length.out=lns[jj]))
                  tix = tix+lns[jj];
                }
              } else if (jj==3){
                if (sum(ts[[3]]!=0)==1 & replRE_1){
                  # If onlye 1 entry in the TS, i.e. checking one particular interval for special rate.
                  # Here first and last are not in, so the
                  if (pfix==1){
                    ltmpfun[[3]] <- eval(substitute(function(x){x[3] +
                        ts[[j1]]*x[alphinx[[j1]]]+ c(0,
                                                     x[seq(stix,length.out=(which(ts[[3]]>0)-2))],0,
                                                     x[seq(stix+which(ts[[3]]>0)-2,length.out=(length(dts)-which(ts[[3]]>0)-1))],
                                                     0)},list(j1=jj,stix=tix)))
                    ptmpfun[[3]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-3))],0,exp(x[3+cumsum(RE)[3]]),log=T))}, list(stix=tix)))
                    reix[[jj]] <- c(seq(tix,length.out=lns[jj]-3))
                    tix = tix+lns[jj]-3;
                  } else if (pfix==2){
                    ltmpfun[[3]] <- eval(substitute(function(x){x[3] +
                        ts[[3]]*x[alphinx[[3]]]+
                        c(x[stix],
                          x[seq(stix,length.out=(which(ts[[3]]>0)-2))],0,
                          x[seq(stix+which(ts[[3]]>0)-2,length.out=(length(dts)-which(ts[[3]]>0)-1))],
                          x[stix+length(dts)-4])},list(j1=jj,stix=tix)))
                    ptmpfun[[3]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-3))],0,exp(x[3+cumsum(RE)[3]]),log=T))}, list(stix=tix)))
                    reix[[jj]] <- c(seq(tix,length.out=lns[jj]-3))
                    tix = tix+lns[jj]-3;
                  }
                } else {
                  # HAVE NOT IMPLEMENTED REMOVAL OF RE IF ONLY 1 STAGE WITH EXCP rate. Problem there is that
                  # there might be more than one timeseries with exceptional rates, which makes the stpwise seq above not work.
                  # This should now be done, see above.
                  if (pfix==1){
                    ltmpfun[[3]] <- eval(substitute(function(x){
                      x[3] +
                        ts[[j1]]*x[alphinx[[j1]]]+ c(0,x[seq(stix,length.out=(length(dts)-2))],0)},list(j1=jj,stix=tix)))
                    ptmpfun[[3]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-2))],0,exp(x[3+cumsum(RE)[3]]),log=T))}, list(stix=tix)))
                    reix[[jj]] <- c(seq(tix,length.out=lns[jj]-2))
                    tix = tix+lns[jj]-2;
                  } else if (pfix==2){
                    ltmpfun[[3]] <- eval(substitute(function(x){x[3] +
                        ts[[j1]]*x[alphinx[[j1]]]+ c(x[stix],x[seq(stix,length.out=(length(dts)-2))],x[stix+length(dts)-3])},list(j1=jj,stix=tix)))
                    ptmpfun[[3]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-2))],0,exp(x[3+cumsum(RE)[3]]),log=T))}, list(stix=tix)))
                    reix[[jj]] <- c(seq(tix,length.out=lns[jj]-2))
                    tix = tix+lns[jj]-2;
                  }
                }
              }
            } else {
              # rowSums(sapply(1:dim(ts[[jj]])[1],function(ii){x[alphinx[[jj]][ii]]*ts[ii,]}))
              # more covariates for same parameter.
              if (jj==1){
                ltmpfun[[1]] <-   eval(substitute(function(x){x[1] +
                    rowSums(sapply(1:dim(ts[[j1]])[1],function(ii){x[alphinx[[j1]][ii]]*ts[[j1]][ii,]}))+
                    c(  x[seq(stix,                 length.out=(length(dts)-1))])},list(j1=jj,stix=tix)))
                ptmpfun[[1]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-1))],0,exp(x[3+cumsum(RE)[1]]),log=T))}, list(stix=tix)))
                reix[[jj]] <- c(seq(tix,length.out=lns[jj]))
                tix = tix+lns[jj];
              } else if (jj==2){
                ltmpfun[[2]] <-   eval(substitute(function(x){x[2] + rowSums(sapply(1:dim(ts[[j1]])[1],function(ii){x[alphinx[[j1]][ii]]*ts[[j1]][ii,]}))+ c(  x[seq(stix, length.out=(length(dts)-1))])},list(stix=tix,j1=jj)))
                ptmpfun[[2]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-1))],0,exp(x[3+cumsum(RE)[2]]),log=T))}, list(stix=tix)))
                reix[[jj]] <- c(seq(tix,length.out=lns[jj]))
                tix = tix+lns[jj];
              } else if (jj==3){
                if (pfix==1){
                  ltmpfun[[3]] <-   eval(substitute(function(x){x[3] + rowSums(sapply(1:dim(ts[[j1]])[1],function(ii){x[alphinx[[j1]][ii]]*ts[[j1]][ii,]}))+ c(0,x[seq(stix,length.out=(length(dts)-2))],0)},list(stix=tix,j1=jj)))
                  ptmpfun[[3]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-2))],0,exp(x[3+cumsum(RE)[3]]),log=T))}, list(stix=tix)))
                  reix[[jj]] <- c(seq(tix,length.out=lns[jj]-2))
                  tix = tix+lns[jj]-2;
                } else if (pfix==2){
                  ltmpfun[[3]] <-   eval(substitute(function(x){x[3] + rowSums(sapply(1:dim(ts[[j1]])[1],function(ii){x[alphinx[[j1]][ii]]*ts[[j1]][ii,]}))+ c(x[stix],x[seq(stix,length.out=(length(dts)-2))],x[stix+length(dts)-3])},list(stix=tix,j1=jj)))
                  ptmpfun[[3]] <-   eval(substitute(function(x){sum(dnorm(x[seq(stix,length.out=(length(dts)-2))],0,exp(x[3+cumsum(RE)[3]]),log=T))}, list(stix=tix)))
                  reix[[jj]] <- c(seq(tix,length.out=lns[jj]-2))
                  tix = tix+lns[jj]-2;
                }
              }
            }
          }
        }
      }
    }
    # Can se just augment and add a term to the ltmps?
    n_est <- function(x){ (n[1:length(dts)-1]/
                             (rate2prob(exp(ltmpfun[[3]](x)),dts)[1:(length(dts)-1)]))}
    n_norm <- function(x){( n_est(x)-mean(n_est(x)))/sd(n_est(x))};
    ltmpfun[[4]] = function(x){0}; # only filled with something if DivDep is true.
    ltmpfun[[5]] = function(x){0};
    if (sum(DivDep==T)>0){
      # Diversity dependence in any of the rates

      if (DivDep[1]==T){
        ltmpfun[[4]] = function(x){x[(4+sum(RE))]*n_norm(x)}
      }  else {
        ltmpfun[[4]] = function(x){0};
      }
      if (DivDep[2]==T){
        ltmpfun[[5]] = function(x){x[(4+sum(RE)+DivDep[1]*1)]*n_norm(x)}
      } else {
        ltmpfun[[5]] = function(x){0};
      }
      # ptmpfun[[4]] <- function(x){sum(dnorm(x[seq(4+sum(RE),length.out=sum(DivDep))],0,5,log=T))};
    }

    myll <- function(x){
      # Simple 3 rate model. Estimated on log x. Fecundity
      # is assumed to be poisson-like, i.e. actual speciation/fecundity for an interval
      # of duration dt is fec*dt.
      lfec = ltmpfun[[1]](x) + ltmpfun[[4]](x);#rep(x[1],length(dts)-1);
      lext = ltmpfun[[2]](x) + ltmpfun[[5]](x);
      lsmp = ltmpfun[[3]](x);

      ll <- pradel_unvd_gam(
        ext = sapply(1:(length(dts)-1),function(ii){
          (exp(lext[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2)))-1))/
            (exp(lfec[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))))-exp(lext[ii]))}),
        gam <- sapply(1:(length(dts)-1),function(ii){
          (1-(exp(lext[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2)))-1))/
             (exp(lfec[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))))-exp(lext[ii])))/
            exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))}),
        p   = sapply(1:(length(dts)),function(ii){rate2prob(exp(lsmp[ii]),dts[ii])}),
        u,n,v,d);
      return(ll)
    };

    # Addition 03.01.2019. Implementing an option for defining priors independently for the
    # three global rates. We might want to 'force' some of these rates to be lower...
    # nprs2 <- priorsNorm_Mus;
    # Doesn't help too much, the priors are not very forceful, the likelihood is too strong.

    myf <- function(x){
      # Simple 3 rate model. Estimated on log x. Fecundity
      # is assumed to be poisson-like, i.e. actual speciation/fecundity for an interval
      # of duration dt is fec*dt.
      lfec = ltmpfun[[1]](x) + ltmpfun[[4]](x);#rep(x[1],length(dts)-1);
      lext = ltmpfun[[2]](x) + ltmpfun[[5]](x);
      lsmp = ltmpfun[[3]](x);
      ll <- pradel_unvd_gam(
        ext = sapply(1:(length(dts)-1),function(ii){
          (exp(lext[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2)))-1))/
            (exp(lfec[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))))-exp(lext[ii]))}),
        gam <- sapply(1:(length(dts)-1),function(ii){
          (1-(exp(lext[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2)))-1))/
             (exp(lfec[ii])*((exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))))-exp(lext[ii])))/
            exp((exp(lfec[ii])-exp(lext[ii]))*(dts[ii]/2+dts[ii+1]/2))}),
        p   = sapply(1:(length(dts)),function(ii){rate2prob(exp(lsmp[ii]),dts[ii])}),
        u,n,v,d);
      ll$LogL <- ll$LogL +
        sum(sapply(1:3,function(ii){dnorm(x[ii],nprs2[[ii]][1],nprs2[[ii]][2],log=T)})) +
        sum(unlist(sapply(1:3,function(a){ptmpfun[[a]](x)}))) +
        (sum(RE)>0)*sum(dunif(x[seq(4,length.out=sum(RE))],priun[1],priun[2],log=T)) + #prior for log(sd(RE)), 0 if RE=c(F;F;F)
        (nhps>(4+sum(RE)))*sum(dnorm(x[seq(4+sum(RE),nhps)],nprs1[1],nprs1[2],log=T)); # prior for covariate efSpecTS.

      # SHould the prior for the log(sd(RE)) be dnorm(0,5)
      if (is.infinite(ll$LogL)){ll$LogL = -Inf};
      if (is.na(ll$LogL)){ll$LogL = -Inf};
      if (is.nan(ll$LogL)){ll$LogL = -Inf};
      return(ll$LogL)
      # So which RE is coded in RE=(T,F,T).
      # no of RE's for each group is lns
    }
  }
  out <- list(probfun= myf,likfun=myll,npar = tix-1,ratefunc= ltmpfun,nhps=nhps,aix = alphinx,ts=ts,ts_orig=ts_orig,lns=lns,priorfun=ptmpfun,n_est=n_est,n_norm=n_norm,reix=reix,Obs=Obs,dts=dts,
              date=date(),call=fullcall);
  attr(out,"class")<-"CMR_model"
  return(out)
  # return(list(probfun= myf,likfun=myll,npar = tix-1,ratefunc= ltmpfun,nhps=nhps,aix = alphinx,ts=ts,lns=lns,priorfun=ptmpfun,n_est=n_est,n_norm=n_norm,reix=reix,Obs=Obs,dts=dts))
}
