#' Function to construct two clade CMR-model.
#'
#' This function generates a structure containing the necessary functions for a
#' CMR analysis of two clades. In addition to drivers affecting both clades, their rates of
#' speciation and extinction can depend on the estimated diversity of both clades.
#'
#' Both clades must be observed at each interval and these must be the same for the two clades.
#'
#' @param Obs1 observations clade 1: a matrix with size \emph{number of taxa} by \emph{number of intervals (n)}. Each taxa has a row with 0's (unobserved) and 1's (observed) for each interval in the analysis. Oldest interval is first column.
#' @param Obs2 observations clade 1: a matrix with size \emph{number of taxa} by \emph{number of intervals (n)}. Each taxa has a row with 0's (unobserved) and 1's (observed) for each interval in the analysis. Oldest interval is first column.
#' @param dts a vector of interval durations. Defaults to series of 1's if not given.
#' @param spec1/ext1/samp1 Formulas for clade 1. Estimated diversity of clade 1 is denoted div1 and for clade 2, div2. E.g. spec1 = ~div1+div2+time, will include two impacts on the speciatino rate of clade 1, estimated diversity of clade 1 (div1) and estimated diversity of clade 2 (div2). These can also interact with other drivers, e.g. spec1= ~div1*driver1 + time
#' @param pfix a switch to select which approach is used to solve the identifiability problem in the model. If \emph{pfix = 1}, then the sampling rate in the first and the last intervals are assumed to be equal to the mean sampling rate for the whole period. If \emph{pfix = 2}, then the first two intervals, and the last two intervals have the same sampling rate. Defaults to pfix = 2.
#' @param data data frame of possible external drivers.
#' @return CMR_model
#' @export
#' @examples M1 <- make_BayesCMR_2clades(Obs1,Obs2,dts=rep(2,5))
#'
#'
#'
make_BayesCMR_2clades <- function(Obs1,Obs2,dts=rep(1,dim(Obs1)[2]),
                                 spec1 = ~time, ext1 = ~time, samp1 = ~time,
                                 spec2 = ~time, ext2 = ~time, samp2 = ~time,pfix=2,data=NULL,...){
  # The goal with this is to actually have formulas as inputs.
  # There must be 6 formulas for this model (spec1,mu1,rho1) (spec2,mu2,rho2)
  # Alternatively we could say that if spec2 is not given, its same as spec1?
  # THis function generates a model with two interacting clades.
  # Drivers can also be fed into clades separately (SpecTS_1 v SpecTS_2 etc).
  # DivDepSpecMatrix and DivDepExtMatrix are matrices which define the interactions between clades. If NA, then no interaction. I
  # If speciation rates of clade 1 is impacte by its own and the other clades diversity (with different effects), then
  # the first row should be [1,2], with numbers indicating unique interactions. If clade1 is equally impacted in its speciation rate
  # from the richness of both own and clade 2, then same numbers can be given
  # DivDepSpecMatrix = [1 , NA]
  #           [NA, 1]
  # denotes if interactions are between and within clades on speciation or extinction rates

  # Need to substitute div1/div2 with div OR

  # But how to make this functionality work if we want the same effect of div1 and div2 on spec1? I.e. if spec1 ~ intcp + beta[1]*(div1+div2)?
  # I guess in principle we could do as mcmc_glmm does with bivariate (i.e. explicitly asking if the spec-rates should be different for
  # each clade AND if the drivers also). Well, this wouldn't really work for (div1,div2), but rather inversely that spec1~beta[1] *div2 and spec2~beta[1]~div2
  # I guess we're back at the matrices.

  # OR we could introduce one more 'div' variable? div1, div2 for diversity of clade 1 and 2, and div12 for their sum, i.e. a common parameter.
  # origdata = drivers;
  # data = drivers
  # normalizing data
  ## TO DO: If last entry (stage) is included in the data.frame, must be removed from normalization!!!
  # Split into two data's, one for samp (length(dts)) and one for spec/ext (SE)
  if (!is.null(data)){
    if (dim(data)[2]>1){
      dataSE = as.data.frame(apply(data[-length(dts),] ,2,normfun))
      data   = as.data.frame(apply(data,2,normfun))
    } else {
      dataSE = as.data.frame(normfun(data[-length(dts)]))
      data = as.data.frame(normfun(data))
    }
  } else {
    dataSE = data.frame(tmp=rep(0,length(dts)-1));
    data   = data.frame(tmp=rep(0,length(dts)));# generating bollocks data.frames so it works when there are no drivers, but temporal variability. Used in making the model.matrices below.
  }

  # Perhaps call this w/o any formula for spec/ext, and implement those outside?
  # OR actually feed in data with div1 and div2 (then treated as 'drivers'), and repopulate the data.frames outside?
  # diversity can only drive S/E's, so matrix should loose last entry in drivers.
  # spec1 <- ~ div1*d180_cor + d13C
  # ext1   <- ~time + div2*d13C
  # samp1 <- ~time
  # spec2 <- spec1
  # ext2 <- ext1
  muk1 <- make_unvd(Obs1*1);
  u1 <- muk1$u
  n1 <- muk1$n
  v1 <- muk1$v
  d1 <- muk1$d

  muk2 <- make_unvd(Obs2*1);
  u2 <- muk2$u
  n2 <- muk2$n
  v2 <- muk2$v
  d2 <- muk2$d

  # dts = m1$dts

  data_tmp <- data.frame(cbind(data,div1=rep(0.1,length(dts)),
                               div2=rep(0.2,length(dts)),
                               div12 = rep(0.3,length(dts))))

  m1 <- make_BayesCMR(  Obs1, dts = dts,spec=spec1,ext=ext1,samp=samp1,data=data_tmp,pfix=pfix,...)
  m2 <- make_BayesCMR(  Obs2, dts = dts,spec=spec2,ext=ext2,samp=samp2,data=data_tmp,pfix=pfix,...)
  # the mmspec outputted
  # Data/drivers for S/E
  # data_tmp <- data.frame(cbind(dataSE,div1=rep(0.1,length(dts)-1),
  #
  #                              div2=rep(0.2,length(dts)-1),div12 = rep(0.3,length(dts)-1)))
  inx <- c();
  inx$inx1 <- m1$inx
  inx$inx2 <- m2$inx
  inx$inx2[1:7]  <- sapply(1:length(m2$inx),function(ii){m2$inx[[ii]]+max(unlist(m1$inx))})
  # The ones above are used for rate calcs, the ones below for n_norm etc
  inxs1 <- seq(1,length.out=m1$npar)
  inxs2 <- seq(1 + max(inxs1),length.out=m2$npar)
  # This is the probfun.
  # ltmpfun are m1$rat
  #
  # inx$xinx1 and the stuff that goes into all these are
  # not the same.

  data_tmpxx <- function(x){
    data.frame(dataSE,div1  = m1$n_norm(x[inxs1]),
               div2  = m2$n_norm(x[inxs2]),
               div12 = m1$n_norm(x[inxs1])*m2$n_norm(x[inxs2]))}

  tmpext1  <- update(ext1,~.-time)
  tmpext2  <- update(ext2,~.-time)
  tmpspec1 <- update(spec1,~.-time)
  tmpspec2 <- update(spec2,~.-time)
  tmpsamp1 <- update(samp1,~.-time)
  tmpsamp2 <- update(samp2,~.-time)
  # So the Zs will be the same for samps
  # We need some check on which of these actually have any of the diversities
  # and need to be functions of x.
  # model.matrix(tmpspec1, data_tmpxx(x)) %*% x[inx$inx1$specInx]
  # model.matrix(tmpspec2, data_tmpxx(x)) %*% x[inx$inx2$specInx]
  # model.matrix(tmpext1,data_tmpxx(x)) %*% x[inx$inx1$extInx]
  # model.matrix(tmpext2,data_tmpxx(x)) %*% x[inx$inx2$extInx]

  # making Z matrices for clade 1
  if ("time" %in% (attr(terms(spec1),"term.labels"))){
    # random effects/temporal effect speciation rates
    Zspec1 <- diag(length(dts)-1) # basically the RE for spec
  } else {
    Zspec1 <- NULL
  }
  if ("time" %in% (attr(terms(ext1),"term.labels"))){
    # random effects/temporal effect speciation rates
    Zext1 <- diag(length(dts)-1) # basically the RE for spec
  } else {
    Zext1 <- NULL
  }
  if ("time" %in% (attr(terms(samp1),"term.labels"))){
    # random effects/temporal effect speciation rates
    if (pfix==2){
      # first two and last two sampling rates are linked
      Zsamp1 <- rbind(c(1,rep(0,length(dts)-3)),diag(length(dts)-2),c(rep(0,length(dts)-3),1))
    } else if (pfix==1){
      # first is equal to last is equal to mean, i.e. no RE for these two
      Zsamp1 <- rbind(rep(0,length(dts)-2),diag(length(dts)-2),rep(0,length(dts)-2))
    }
  } else {
    Zsamp1 <- NULL
  }
  # Making Z matrices for clade 2
  if ("time" %in% (attr(terms(spec2),"term.labels"))){
    # random effects/temporal effect speciation rates
    Zspec2 <- diag(length(dts)-1) # basically the RE for spec
  } else {
    Zspec2 <- NULL
  }
  if ("time" %in% (attr(terms(ext2),"term.labels"))){
    # random effects/temporal effect speciation rates
    Zext2 <- diag(length(dts)-1) # basically the RE for spec
  } else {
    Zext2 <- NULL
  }
  if ("time" %in% (attr(terms(samp2),"term.labels"))){
    # random effects/temporal effect speciation rates
    # note, same pfix for both clades)
    if (pfix==2){
      # first two and last two sampling rates are linked
      Zsamp2 <- rbind(c(1,rep(0,length(dts)-3)),diag(length(dts)-2),c(rep(0,length(dts)-3),1))
    } else if (pfix==1){
      # first is equal to last is equal to mean, i.e. no RE for these two
      Zsamp2 <- rbind(rep(0,length(dts)-2),diag(length(dts)-2),rep(0,length(dts)-2))
    }
  } else {
    Zsamp2 <- NULL
  }


  # ("div" %in% substr((attr(terms(spec1),"term.labels")),1,3))
  # Some diversity in this
  # So if no diversity, just use dataSE, if diversity ,use data_tmp(x) as input to model matrix
  # model.matrix(tmpspec1, data_tmpxx(x)) %*% x[inx$inx1$specInx]
  # model.matrix(tmpspec2, data_tmpxx(x)) %*% x[inx$inx2$specInx]
  # model.matrix(tmpext1,data_tmpxx(x)) %*% x[inx$inx1$extInx]
  # model.matrix(tmpext2,data_tmpxx(x)) %*% x[inx$inx2$extInx]


  #### ==== ####
  # defining rate functions
  # Clade1 speciation function
  if (is.null(Zspec1)){
    # no random effects
    if (length(inx$inx1$specInx)==1){
      # if only intercept
      lspecfun1 <- function(x){rep(x[inx$inx1$specInx],length(dts)-1)}
    } else {
      # here be drivers! But no random effects
      if ("div" %in% substr((attr(terms(spec1),"term.labels")),1,3)){
        # diversity part of drivers, need to use function for dataframe
        lspecfun1 <- function(x){
          model.matrix(tmpspec1,data_tmpxx(x)) %*% x[inx$inx1$specInx]}
      } else {
        # if no diversity, use datasE
        lspecfun1 <- function(x){
          model.matrix(tmpspec1,dataSE) %*% x[inx$inx1$specInx]}
      }
    }
  } else {
    # also with random effects
    if (length(inx$inx1$specInx)==1){
      # if only intercept
      lspecfun1 <- function(x){
        rowSums(cbind(model.matrix(tmpspec1,dataSE) %*% x[inx$inx1$specInx], Zspec1 %*% x[inx$inx1$specReInx]))}
    } else {
      # here be drivers! But with random effects
      if ("div" %in% substr((attr(terms(spec1),"term.labels")),1,3)){
        # diversity part of drivers, need to use function for dataframe
        lspecfun1 <- function(x){
          rowSums(cbind(model.matrix(tmpspec1,data_tmpxx(x)) %*% x[inx$inx1$specInx],Zspec1 %*% x[inx$inx1$specReInx]))}
      } else {
        # if no diversity, use datasE
        lspecfun1 <- function(x){
          rowSums(cbind(model.matrix(tmpspec1,dataSE) %*% x[inx$inx1$specInx],Zspec1 %*% x[inx$inx1$specReInx]))}
      }
    }
  }


  # Clade2 speciation function
  if (is.null(Zspec2)){# no random effects
    if (length(inx$inx2$specInx)==1){# if only intercept
        lspecfun2 <- function(x){rep(x[inx$inx2$specInx],length(dts)-1)}
    } else {# here be drivers! But no random effects
      if ("div" %in% substr((attr(terms(spec2),"term.labels")),1,3)){# diversity part of drivers, need to use function for dataframe
        lspecfun2 <- function(x){model.matrix(tmpspec2,data_tmpxx(x)) %*% x[inx$inx2$specInx]}
      } else {# if no diversity, use datasE
        lspecfun2 <- function(x){model.matrix(tmpspec2,dataSE) %*% x[inx$inx2$specInx]}
      }
    }
  } else { # also with random effects
    if (length(inx$inx2$specInx)==1){ # if only intercept
        lspecfun2 <- function(x){ rowSums(cbind(model.matrix(tmpspec2,dataSE)       %*% x[inx$inx2$specInx], Zspec2 %*% x[inx$inx2$specReInx]))}
    } else { # here be drivers! But WITH random effects
      if ("div" %in% substr((attr(terms(spec2),"term.labels")),1,3)){ # diversity part of drivers, need to use function for dataframe
        lspecfun2 <- function(x){rowSums(cbind(model.matrix(tmpspec2,data_tmpxx(x)) %*% x[inx$inx2$specInx], Zspec2 %*% x[inx$inx2$specReInx]))}
      } else {# if no diversity, use datasE
        lspecfun2 <- function(x){rowSums(cbind(model.matrix(tmpspec2,dataSE)        %*% x[inx$inx2$specInx], Zspec2 %*% x[inx$inx2$specReInx]))}
      }
    }
  }


  # Clade1 extinction function
  if (is.null(Zext1)){# no random effects
    if (length(inx$inx1$extInx)==1){# if only intercept
      lextfun1 <- function(x){rep(x[inx$inx1$extInx],length(dts)-1)}
    } else {# here be drivers! But no random effects
      if ("div" %in% substr((attr(terms(ext1),"term.labels")),1,3)){# diversity part of drivers, need to use function for dataframe
        lextfun1 <- function(x){model.matrix(tmpext1,data_tmpxx(x)) %*% x[inx$inx1$extInx]}
      } else {# if no diversity, use datasE
        lextfun1 <- function(x){model.matrix(tmpext1,dataSE) %*% x[inx$inx1$extInx]}
      }
    }
  } else { # also with random effects
    if (length(inx$inx1$extInx)==1){ # if only intercept
      lextfun1 <- function(x){ rowSums(cbind(model.matrix(tmpext1,dataSE)       %*% x[inx$inx1$extInx], Zext1 %*% x[inx$inx1$extReInx]))}
    } else { # here be drivers! But WITH random effects
      if ("div" %in% substr((attr(terms(ext1),"term.labels")),1,3)){ # diversity part of drivers, need to use function for dataframe
        lextfun1 <- function(x){rowSums(cbind(model.matrix(tmpext1,data_tmpxx(x)) %*% x[inx$inx1$extInx], Zext1 %*% x[inx$inx1$extReInx]))}
      } else {# if no diversity, use datasE
        lextfun1 <- function(x){rowSums(cbind(model.matrix(tmpext1,dataSE)        %*% x[inx$inx1$extInx], Zext1 %*% x[inx$inx1$extReInx]))}
      }
    }
  }

  # Clade2 extinction function
  if (is.null(Zext2)){# no random effects
    if (length(inx$inx2$extInx)==1){# if only intercept
      lextfun2 <- function(x){rep(x[inx$inx2$extInx],length(dts)-1)}
    } else {# here be drivers! But no random effects
      if ("div" %in% substr((attr(terms(ext2),"term.labels")),1,3)){# diversity part of drivers, need to use function for dataframe
        lextfun2 <- function(x){model.matrix(tmpext2,data_tmpxx(x)) %*% x[inx$inx2$extInx]}
      } else {# if no diversity, use datasE
        lextfun2 <- function(x){model.matrix(tmpext2,dataSE) %*% x[inx$inx2$extInx]}
      }
    }
  } else { # also with random effects
    if (length(inx$inx2$extInx)==1){ # if only intercept
      lextfun2 <- function(x){ rowSums(cbind(model.matrix(tmpext2,dataSE)       %*% x[inx$inx2$extInx], Zext2 %*% x[inx$inx2$extReInx]))}
    } else { # here be drivers! But WITH random effects
      if ("div" %in% substr((attr(terms(ext2),"term.labels")),1,3)){ # diversity part of drivers, need to use function for dataframe
        lextfun2 <- function(x){rowSums(cbind(model.matrix(tmpext2,data_tmpxx(x)) %*% x[inx$inx2$extInx], Zext2 %*% x[inx$inx2$extReInx]))}
      } else {# if no diversity, use datasE
        lextfun2 <- function(x){rowSums(cbind(model.matrix(tmpext2,dataSE)        %*% x[inx$inx2$extInx], Zext2 %*% x[inx$inx2$extReInx]))}
      }
    }
  }

  # Clade1 sampling function
  if (is.null(Zsamp1)){# no random effects
    if (length(inx$inx1$sampInx)==1){# if only intercept
        lsampfun1 <- function(x){rep(x[inx$inx1$sampInx],length(dts))}
    } else {# here be drivers! But no random effects and diversity cannot be one .Use data, dataSE, since samps have
      # BUT, with drivers and SE's, we are assuming that the RE's are identical for the last /first two stages. But not that the rates are identical!
        lsampfun1 <- function(x){model.matrix(tmpsamp1,data) %*% x[inx$inx1$sampInx]}
    }
  } else { # also with random effects
    if (length(inx$inx1$sampInx)==1){ # if only intercept
        lsampfun1 <- function(x){ rowSums(cbind(model.matrix(tmpsamp1,data)      %*% x[inx$inx1$sampInx], Zsamp1 %*% x[inx$inx1$sampReInx]))}
    } else { # here be drivers! But WITH random effects
        lsampfun1 <- function(x){ rowSums(cbind(model.matrix(tmpsamp1,data)        %*% x[inx$inx1$sampInx], Zsamp1 %*% x[inx$inx1$sampReInx]))}
    }
  }

  # Clade2 sampling function
  if (is.null(Zsamp2)){# no random effects
    if (length(inx$inx2$sampInx)==1){# if only intercept
      lsampfun2 <- function(x){rep(x[inx$inx2$sampInx],length(dts))}
    } else {# here be drivers! But no random effects and diversity cannot be one .Use data, dataSE, since samps have
      # BUT, with drivers and SE's, we are assuming that the RE's are identical for the last /first two stages. But not that the rates are identical!
      lsampfun2 <- function(x){model.matrix(tmpsamp2,data) %*% x[inx$inx2$sampInx]}
    }
  } else { # also with random effects
    if (length(inx$inx2$sampInx)==1){ # if only intercept
      lsampfun2 <- function(x){ rowSums(cbind(model.matrix(tmpsamp2,data)      %*% x[inx$inx2$sampInx], Zsamp2 %*% x[inx$inx2$sampReInx]))}
    } else { # here be drivers! But WITH random effects
      lsampfun2 <- function(x){ rowSums(cbind(model.matrix(tmpsamp2,data)        %*% x[inx$inx2$sampInx], Zsamp2 %*% x[inx$inx2$sampReInx]))}
    }
  }

  #### ==== ####

# Making the probability density
  myll <- function(x){
    lfec1 = lspecfun1(x);
    lext1 = lextfun1(x);
    lsmp1  = lsampfun1(x);
    lfec2  = lspecfun2(x);
    lext2  = lextfun2(x);
    lsmp2  = lsampfun2(x);

    ll1 <- pradel_unvd_gam(
      ext = sapply(1:(length(dts)-1),function(ii){
          (exp(lext1[ii])*((exp((exp(lfec1[ii])-exp(lext1[ii]))*dts[ii]))-1))/
          (exp(lfec1[ii])*((exp((exp(lfec1[ii])-exp(lext1[ii]))*dts[ii])))-exp(lext1[ii]))}),
      gam <- sapply(1:(length(dts)-1),function(ii){
        (1-(exp(lext1[ii])*((exp((exp(lfec1[ii])-exp(lext1[ii]))*dts[ii]))-1))/
           (exp(lfec1[ii])*((exp((exp(lfec1[ii])-exp(lext1[ii]))*dts[ii])))-exp(lext1[ii])))/
          exp((exp(lfec1[ii])-exp(lext1[ii]))*dts[ii])}),
      p   = sapply(1:(length(dts)),function(ii){rate2prob(exp(lsmp1[ii]),dts[ii])}),
      u1,n1,v1,d1);

    ll2 <- pradel_unvd_gam(
      ext = sapply(1:(length(dts)-1),function(ii){
          (exp(lext2[ii])*((exp((exp(lfec2[ii])-exp(lext2[ii]))*dts[ii]))-1))/
          (exp(lfec2[ii])*((exp((exp(lfec2[ii])-exp(lext2[ii]))*dts[ii])))-exp(lext2[ii]))}),
      gam <- sapply(1:(length(dts)-1),function(ii){
        (1-(exp(lext2[ii])*((exp((exp(lfec2[ii])-exp(lext2[ii]))*dts[ii]))-1))/
           (exp(lfec2[ii])*((exp((exp(lfec2[ii])-exp(lext2[ii]))*dts[ii])))-exp(lext2[ii])))/
          exp((exp(lfec2[ii])-exp(lext2[ii]))*dts[ii])}),
      p   = sapply(1:(length(dts)),function(ii){rate2prob(exp(lsmp2[ii]),dts[ii])}),
      u2,n2,v2,d2);

    if (is.infinite(ll1$LogL)){ll1$LogL = -Inf};
    if (is.na(ll1$LogL)){ll1$LogL = -Inf};
    if (is.nan(ll1$LogL)){ll1$LogL = -Inf};
    if (is.infinite(ll2$LogL)){ll2$LogL = -Inf};
    if (is.na(      ll2$LogL)){ll2$LogL = -Inf};
    if (is.nan(     ll2$LogL)){ll2$LogL = -Inf};
    return(list(ll1=ll1,ll2=ll2))
  }

  probfun <- function(x){myll(x)$ll1$LogL + myll(x)$ll2$LogL +
      sum(unlist(m1$priorf(x[inxs1])))+
      sum(unlist(m2$priorf(x[inxs2])))}
  # THis is strictly speaking not true. If div12 is in both
  # Actually this points to a problem; interactions can not be the same for the two clades;
  # i.e. spec1 ~ b*div1 and spec2 ~b*div2, will not work. spec1~b*div12 and spec2~div12 will, however.
  # THis is actually the same for drivers; they can not be assumed to have
  # the same impact on both clades. This could potentially be done by
  # some meedling of the index array, but leave for later.
  out <- list(probfun = probfun, likfun = myll,Clade1Mod = m1, Clade2Mod = m2,
              lspecfun1 = lspecfun1,lextfun1=lextfun1,lsampfun1=lsampfun1,
              lspecfun2 = lspecfun2,lextfun2=lextfun2,lsampfun2=lsampfun2,
              inx = inx,npar = max(unlist(inx)),date=date());

  attr(out,"class")<-"CMR_model"
  return(out)
}


