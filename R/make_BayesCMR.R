#' Function to construct one clade CMR-model.
#'
#' This function generates a structure containing the necessary functions for a CMR analysis of a fossil dataset.
#'
#' @param Obs a matrix with size \emph{number of taxa} by \emph{number of intervals (n)}. Each taxa has a row with 0's (unobserved) and 1's (observed) for each interval in the analysis. Oldest interval is first column.
#' @param dts a vector of interval durations. Defaults to a series of 1's if not supplied. Oldest interval is first entry.
#' @param spec/ext/samp are formulas specifying the particular model. Use 'time' for temporally variable rates, 'div' for diversity dependence. Other drivers should have names as entries into \emph{data}
#' @param data a data frame with putative drivers.
#' @param pfix a switch to select which approach is used to solve the identifiability problem in the model. If \emph{pfix = 1}, then the sampling rate in the first and the last intervals are assumed to be equal to the mean sampling rate for the whole period. If \emph{pfix = 2}, then the first two intervals, and the last two intervals have the same sampling rate. Defaults to pfix = 2.
#' @param priorsNorm_Cov sets the parameters for the normal prior for the impact of the drivers. Defaults to Norm(mu=0,sd=2).
#' @param priorsNorm_Mus sets the parameters for the normal prior for the (log) mean speciation, extinction and sampling rates. Given as a list of three entries with (mu,sd). Defaults to rep(list(c(-4,4)),3), i.e. all normal priors with mu=-4 and sd=3.
#' @param priorsUnif sets the parameters for the uniform prior for the log variances of the random efSpecTS. Same for all three variances, dunif(-3,3).
#' @param modeltype sets the type parameterization of the rates of speciation and extinction. 'I' and 'II' specifies that the probabilities of extinction and seniority are calculated based on actual durations of intervals, whereas 'III' and 'IV' uses the differences between mid-points of intervals. 'I' and 'III' uses a basic implementation of probability of extinction as 1 - exp(-rate * dt). Models 'II' and 'IV' use a more detailed functional relationship between rate and probability, following Raup (1984):  (µ*((exp((lambda-µ)*dt))-1))/(lambda*((exp((lambda-µ)*dt)))-µ). Modeltype 'V' implements a temporally uninformed model, estimating the probabilities directly (or rather their logit).
#' @return The function returns an object of class CMR_model. This most important output is out$probfun which is a function for the posterior. This output can be fed directly into the sampler \link{MCMC_CMR}. It also returns \emph{Obs} and \emph{dts} as well as the full call (exlcluding default settings), \emph{call}
#' @export
#' @examples mod1 <- make.BayesCMR(Obs)
#' fit <- MCMC_CMR(mod1)
#' matplot(fit$Chain[,1:3],type="l")

make_BayesCMR <- function(Obs,dts=rep(1,dim(Obs)[2]),
                              spec =  ~ 1,
                              ext  =  ~ 1,
                              samp =  ~ 1,
                              data = NULL,
                              pfix = 2,
                              priorMu  = dnorm,
                              priorCov = dnorm,
                              priorStd = dunif,
                              priorPars = list(list(-4,4),list(0,2),list(-3,3)),
                              replRE_1 = F,modeltype='IV'){

  # Model generating function for a Compadre analysis.
  # Minimum input is a matrix of observed/unobserved of dimensions
  # taxa by temporal interval. dts is vector of interval
  # durations (default sets to 1).

  # Dts are in this version difference between mid-poitns of intervals for spec/ext and true dts for samp.
  dts_se <- dts[-1]/2 + dts[-length(dts)]/2

  # Current version implements paraclade extinction and growth from Raup 1984
  # used to get seniority and extinction. Current version 26.03.2019

  ## INPUT CHECKS
  fullcall <- deparse(match.call());
  if (length(dts) != dim(Obs)[2]){
    stop("Number of intervals in *dts* does not match dimension of Obs")
  }
  if (!is.null(data)){
    # If drivers are given
    if (NROW(data) != length(dts)){
      stop("Size of data frame with drivers does not match dimension of Obs.")
    }
  }

  undv <- make_unvd(Obs);
  u = undv$u;
  v = undv$v;
  n = undv$n;
  d = undv$d;

  prior <- list(); #making all prior terms part of this list.
  prix = 1; # ticker for prior parts.

  origdata = data;
  # normalizing data
  # If last entry (stage) is included in the data.frame, must be removed from normalization for the Speciation/extinction drivers
  # Split into two data's, one for samp (length(dts)) and one for spec/ext (SE)
  if (!is.null(data)){
    if (dim(data)[2]>1){
      dataSE = as.data.frame(apply(data[-length(dts),] ,2,normfun))
      data   = as.data.frame(apply(data,2,normfun))
    } else {
      dataSE = as.data.frame(normfun(data[-length(dts)]))
      dataSE = as.data.frame(normfun(data))
    }
  } else {
    dataSE = data.frame(tmp=rep(0,length(dts)-1));
    data   = data.frame(tmp=rep(0,length(dts)));
    # generating bollocks data.frames so it works when there are no drivers, but temporal variability. Used in making the model.matrices below.
  }
  updmmx <- c(F,F,F) # if updating the model matrices for any X (i.e. if diversity is a driver)
  # If mmxxx needs to be updated, i.e. if diversity is a term in any of the rates.
  # SO this fails if not data given and div-dep, since data then doesn't exist properly.
  if ("div" %in% (attr(terms(spec),"term.labels"))){
    if (is.null(dataSE)){
      dataSE = data.frame(div=rep(0,length(dts)-1))
    } else {
      dataSE$div = rep(0,dim(dataSE)[1])
    }
    updmmx[1] <- TRUE
    # Perhaps have these as lists with indexes instead? If there is an interaction, several rows (not just $div) must be updated!
  }
  if ("div" %in% (attr(terms(ext),"term.labels"))){
    if (is.null(dataSE)){
      dataSE = data.frame(div=rep(0,length(dts)-1))
    } else {
      dataSE$div = rep(0,dim(dataSE)[1])
    }
    updmmx[2] <- TRUE
  }


  if ("div" %in% (attr(terms(samp),"term.labels"))){
    # This is not possible... we can not have diversity estimated affecting the sampling of 'diversity' itself. Throw error
    error('Model with diversity impact on sampling is not possible.')
  }
  # Making design matrixes for rates. THESE ARE TEMPORARY if diversity is a term. If not they are unchanging.
  tmpspec <- update(spec,~.-time) # taking out TIME, since time is the Random effects.
  mmspec  <- model.matrix(tmpspec,dataSE)
  mmspec_f <- function(x){
    tmp <- dataSE;
    tmp$div <- n_norm(x)
    return(model.matrix(tmpspec,tmp))}
  tmpext  <- update(ext,~.-time) # taking out TIME, since time is the Random effects.
  mmext   <- model.matrix(tmpext,dataSE)
  mmext_f <- function(x){
    tmp <- dataSE;
    tmp$div <- n_norm(x)
    return(model.matrix(tmpext,tmp))}
  tmpsamp <- update(samp,~.-time) # taking out TIME, since time is the Random effects.
  mmsamp  <- model.matrix(tmpsamp,data)


  RE = c(F,F,F)
  if ("time" %in% (attr(terms(spec),"term.labels"))){
    # random effects/temporal effect speciation rates
    Zspec <- diag(length(dts)-1) # basically the RE for spec
    RE[1] = TRUE
  } else {
    Zspec <- NULL
  }

  if ("time" %in% (attr(terms(ext),"term.labels"))){
    # random effects/temporal effect speciation rates
    Zext <- diag(length(dts)-1) # basically the RE for spec
    RE[2] = TRUE
  } else {
    Zext <- NULL
  }

  if ("time" %in% (attr(terms(samp),"term.labels"))){
    # random effects/temporal effect speciation rates
    # Making 'RE' matrices.
    RE[3] = TRUE
    if (pfix==2){
      # first two and last two sampling rates are linked
      Zsamp <- rbind(c(1,rep(0,length(dts)-3)),diag(length(dts)-2),c(rep(0,length(dts)-3),1))
    } else if (pfix==1){
      # first is equal to last is equal to mean, i.e. no RE for these two
      Zsamp <- rbind(rep(0,length(dts)-2),diag(length(dts)-2),rep(0,length(dts)-2))
    }
  } else {
    Zsamp <- NULL
  }

  # cbind(mmsamp*c(0.3), Zsamp*runif(length(dts)))
  # THIS MIGHT BE DIFFERENT FRMO BEFORE, NOW DRIVERS
  # are BEfore variances in the x-array.
  alphinx= list();
  alphinx$specInx   = c(1,seq(4,        length.out=(dim(mmspec)[2]-1)))
  alphinx$extInx    = c(2,seq(max(c(3,alphinx$specInx))+1,length.out=(dim(mmext)[2]-1)))
  alphinx$sampInx   = c(3,seq(max(c(3,alphinx$specInx,alphinx$extInx))+1,length.out=(dim(mmsamp)[2]-1)))


  alphinx$varInx    = seq(max(unlist(alphinx)+1),length.out=sum(RE))
  alphinx$specReInx = seq(max(unlist(alphinx)+1),length.out=RE[1]*(length(dts)-1))
  alphinx$extReInx  = seq(max(unlist(alphinx)+1),length.out=RE[2]*(length(dts)-1))
  alphinx$sampReInx = seq(max(unlist(alphinx)+1),length.out=RE[3]*(length(dts)-2))

  names(alphinx$specInx) <- colnames(mmspec)
  names(alphinx$extInx)  <- colnames(mmext)
  names(alphinx$sampInx) <- colnames(mmsamp)

  npar <- dim(mmspec)[2] + dim(mmext)[2] + dim(mmsamp)[2] + sum(RE*1) + RE[1]*(length(dts)-1) + RE[2]*(length(dts)-1) + RE[3]*(length(dts)-2)
  # This is number of parameters in total.
  # sum(DivDep*1) + nSmpTS+nExtTS+nSpecTS + sum(Driv_x_Div_Spec) + sum(Driv_x_Div_Ext); # number of mean+hyperparameters
  # nTS  <- c(nSpecTS,nExtTS,nSmpTS);
  # ALL x-arrays endd with the RE's, so all indexes below here that are not the first 3+sum(RE) sohuld have sum(DivDep) added to them.
  # All these have priors dnorm(0,5), the Sd[RE] are log(sd) really.
  # [�_fec,�_ext,�samp] &
  # (possibly)[st_RE_fec,st_RE_ext,st_RE_samp] &
  # possibly [alphaTS1,alphaTS2 etc]
  #
  #
  if (is.null(Zsamp)){
    # If no RE for sampling
    if (dim(mmsamp)[2]==1){
      # if only intercept
      lsampfun <- function(x){rep(x[3],length(dts))}
    } else {
      # here with drivers
      lsampfun <- function(x){rowSums(cbind(mmsamp %*% x[alphinx$sampInx]))}
    }
  } else {
    # if RE's
    lsampfun <- function(x){rowSums(cbind(mmsamp %*% x[alphinx$sampInx], Zsamp %*% x[alphinx$sampReInx]))}
  }

  # making div function
  n_est <- function(x){ (n[1:length(dts)-1]/
                           (rate2prob(exp(lsampfun(x)),dts)[1:(length(dts)-1)]))}
  n_norm <- function(x){( n_est(x)-mean(n_est(x)))/sd(n_est(x))};


  if (is.null(Zspec)){# If no RE for speciation
    if (dim(mmspec)[2]==1){# if only intercept
      lspecfun <- function(x){rep(x[1],length(dts)-1)}
    } else {        # here with drivers
      if (updmmx[1]){
        lspecfun <- function(x){rowSums(cbind(mmspec_f(x) %*% x[alphinx$specInx]))       }
      } else {
        lspecfun <- function(x){rowSums(cbind(mmspec %*% x[alphinx$specInx]))}
      }
    }
  } else {# if RE's
    if (updmmx[1]){
      lspecfun <- function(x){
        return(rowSums(cbind(mmspec_f(x) %*% x[alphinx$specInx], Zspec %*% x[alphinx$specReInx])))        }
    } else {
      lspecfun <- function(x){rowSums(cbind(mmspec %*% x[alphinx$specInx], Zspec %*% x[alphinx$specReInx]))}
    }
  }


  if (is.null(Zext)){# If no RE for extinction
    if (dim(mmext)[2]==1){# if only intercept
      lextfun <- function(x){rep(x[2],length(dts)-1)}
    } else {        # here with drivers
      if (updmmx[2]){
        # If div-dep, fill in with divers
        lextfun <- function(x){rowSums(cbind(mmext_f(x) %*% x[alphinx$extInx]))       }
      } else {
        lextfun <- function(x){rowSums(cbind(mmext_f(x) %*% x[alphinx$extInx]))}
      }
    }
  } else {# if RE's
    # Should also be here with drivers, and divdep interactions and all
    # lextfun <- function(x){rowSums(cbind(mmext %*% x[alphinx$extInx], Zext %*% x[alphinx$extReInx]))}
    if (updmmx[2]){
      # If div-dep, fill in with divers
      lextfun <- function(x){rowSums(cbind(mmext_f(x) %*% x[alphinx$extInx], Zext %*% x[alphinx$extReInx]))}
    } else {
      lextfun <- function(x){rowSums(cbind(mmext %*% x[alphinx$extInx], Zext %*% x[alphinx$extReInx]))}
    }
  }

  # Defining the likelihoodfunction.
  myll <- function(x){
    lfec = lspecfun(x);
    lext = lextfun(x);
    lsmp = lsampfun(x);
    if (modeltype=='I'){
      extin <- 1-exp(-exp(lext)*dts[-length(dts)])
      gamin <- (exp(-exp(lext)*dts[-length(dts)]))/(exp((exp(lfec)-exp(lext))*dts[-length(dts)]))
      sampin <- rate2prob(exp(lsmp),dts)
    } else if (modeltype == 'II'){
      extin <-    (exp(lext)*((exp((exp(lfec)-exp(lext))*dts[-length(dts)]))-1))/        (exp(lfec)*((exp((exp(lfec)-exp(lext))*dts[-length(dts)])))-exp(lext))
      gamin <- (1-(exp(lext)*((exp((exp(lfec)-exp(lext))*dts[-length(dts)]))-1))/        (exp(lfec)*((exp((exp(lfec)-exp(lext))*dts[-length(dts)])))-exp(lext)))/        exp((exp(lfec)-exp(lext))*dts[-length(dts)])
      sampin <- rate2prob(exp(lsmp),dts)
    } else if (modeltype == 'III'){
      extin <- 1-exp(-exp(lext)*dts_se)
      gamin <- (exp(-exp(lext)*dts_se))/(exp((exp(lfec)-exp(lext))*dts_se))
      sampin <- rate2prob(exp(lsmp),dts)
    } else if (modeltype == 'IV'){
      extin <- (exp(lext)*((exp((exp(lfec)-exp(lext))*dts_se))-1))/        (exp(lfec)*((exp((exp(lfec)-exp(lext))*dts_se)))-exp(lext))
      gamin <- (1-(exp(lext)*((exp((exp(lfec)-exp(lext))*dts_se))-1))/        (exp(lfec)*((exp((exp(lfec)-exp(lext))*dts_se)))-exp(lext)))/        exp((exp(lfec)-exp(lext))*dts_se)
      sampin <- rate2prob(exp(lsmp),dts)
    } else if (modeltype == 'V'){
      extin <- exp(lext)/(exp(lext)+1)
      gamin <- exp(lfec)/(exp(lfec)+1)
      sampin <- exp(lsmp)/(exp(lsmp)+1)
    }
    ll <- pradel_unvd_gam(
      ext = extin,
      gam = gamin,
      p   = sampin,
      u,n,v,d);
    if (is.infinite(ll$LogL)){ll$LogL = -Inf};
    if (is.na(ll$LogL)){ll$LogL = -Inf};
    if (is.nan(ll$LogL)){ll$LogL = -Inf};
    return(ll)

  }
  inxDriv = seq(4,length.out=dim(mmspec)[2]+ dim(mmext)[2]+dim(mmsamp)[2]-3) # index into 'drivers' non intercepts.
  if (any(RE)){

    prfun <- function(x){unlist(c(sapply(1:3,function(ii){do.call(priorMu, c(x[ii],priorPars[[1]],log=T))}),
                                  sapply(inxDriv,function(ii){do.call(priorCov,c(x[ii],priorPars[[2]],log=T))}),
                                  sapply(alphinx$varInx,function(ii){do.call(priorStd,c(x[ii],priorPars[[3]],log=T))}),
                                  sapply(1:sum(RE),function(ii){dnorm(   x[alphinx[[4+which(RE)[ii]]]],0,exp(x[alphinx$varInx[ii]]),log=T)})))}
  } else { # no random effects, i.e. no varinx
    prfun <- function(x){unlist(c(sapply(1:3,function(ii){do.call(priorMu, c(x[ii],priorPars[[1]],log=T))}),
                                  sapply(inxDriv,function(ii){do.call(priorCov,c(x[ii],priorPars[[2]],log=T))})))}
  }

  probfun <- function(x){myll(x)$LogL+sum(unlist(prfun(x)))}

  out <- list(probfun= probfun,likfun=myll,priorf  = prfun,
              npar = npar,test = 'dummy',
              specfun = lspecfun,extfun=lextfun,sampfun=lsampfun,inx = alphinx,driverinx = inxDriv,
              n_est=n_est,n_norm=n_norm,Obs=Obs,dts=dts,
              date=date(),call=fullcall,modeltype=modeltype,
              mmsamp=mmsamp,mmspec=mmspec_f,mmext=mmext_f,spec = spec,ext = ext,samp=samp)
  # Zsmp = Zsamp, Zspc = Zspec, Zext = Zext,dataSE = dataSE,
  # inxDriv = inxDriv,RE = RE);
  attr(out,"class")<-"CMR_model"
  return(out)


}
