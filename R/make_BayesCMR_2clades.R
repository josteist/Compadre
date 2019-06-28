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
#' @param DivDepSpecMatrix matrix of interactions between clades on speciation rates. Given as unique integers, starting at 1. For only diversity dependence of speciation rates within clades = matrix(c(1,NA,NA,2),nrow=2). For only between clade interactions = matrix(c(NA,1,2,NA),nrow=2). Similarly for DivDepExtMatrix, starting at the ma(DivDepSpecMatrix)+1
#' @param SpecTS_1/2 drivers for speciation rates for clades 1/2
#' @param ExtTS_1/2 drivers for extinction rates for clades 1/2
#' @param SmpTS_1/2 drivers for sampling rates
#' @param Driv_x_Div_Spec_1/2 whether or not speciation drivers interact with diversity in clade 1/2
#' @param Driv_x_Div_Spec_1/2 whether or not extinction drivers interact with diversity in clade1/2. TRUE or FALSE of length equal to number of drivers.
#' @return this function yields whatnow?
#' @export
#' @examples M1 <- make_BayesCMR_2clades(Obs1,Obs2,dts=rep(2,5))
#'
#'
#'
make_BayesCMR_2clades <- function(Obs1,Obs2,dts=rep(1,dim(Obs1)[2]),
                                  RE=c(T,T,T),
                                  DivDepSpecMatrix=matrix(rep(NA,4),nrow=2),
                                  DivDepExtMatrix=matrix(rep(NA,4),nrow=2),
                                  SpecTS_1=NULL,ExtTS_1=NULL,SmpTS_1=NULL,
                                  SpecTS_2=NULL,ExtTS_2=NULL,SmpTS_2=NULL,
                                  Driv_x_Div_Spec_1=NULL,Driv_x_Div_Ext_1=NULL,
                                  Driv_x_Div_Spec_2=NULL,Driv_x_Div_Ext_2=NULL,...){
  # THis function generates a model with two interacting clades.
  # Drivers can also be fed into clades separately (SpecTS_1 v SpecTS_2 etc).
  # DivDepSpecMatrix and DivDepExtMatrix are matrices which define the interactions between clades. If NA, then no interaction. I
  # If speciation rates of clade 1 is impacte by its own and the other clades diversity (with different effects), then
  # the first row should be [1,2], with numbers indicating unique interactions. If clade1 is equally impacted in its speciation rate
  # from the richness of both own and clade 2, then same numbers can be given
  # DivDepSpecMatrix = [1 , NA]
  #           [NA, 1]
  # denotes if interactions are between and within clades on speciation or extinction rates
  m1 <- make_BayesCMR(  Obs1, dts = dts,RE=RE,SpecTS = SpecTS_1,ExtTS= ExtTS_1,SmpTS = SmpTS_1,Driv_x_Div_Spec=Driv_x_Div_Spec_1,Driv_x_Div_Ext=Driv_x_Div_Ext_1,...)
  m2 <- make_BayesCMR(  Obs2, dts = dts,RE=RE,SpecTS = SpecTS_2,ExtTS= ExtTS_2,SmpTS = SmpTS_2,Driv_x_Div_Spec=Driv_x_Div_Spec_2,Driv_x_Div_Ext=Driv_x_Div_Ext_2,...)

  # This is the probfun.
  # ltmpfun are m1$rat
  #
  # TO DO (short term)
  # - make it allow for drivers. DONE
  # - expand diversity effect switches to be matrices for Spec and Ext DONE

  # Possible bugs. It seems like there is a higher risk of non-applicable initial values for this model. Initialize narrowly around 0 and it seems to work.
  # Problem is when some of the rates blow up, which is more likely the more 'impact's with random initial values for parameters affecting the rates.

  # TO DO (long term)
  # - allow for non-overlapping temporal spans; i.e. that one clade has 0 obs and no pars in some intervals
  # So lets just set the X's to be after eachother, first ALL the
  # x's for m1, then m2, and lastly the 1-4 pars for their interactions.
  # This will ONLY work for totally overlapping intervals of observation.
  # NOTE THAT THE MAKE_BAYESCALL was W/O diversity dependence internally, but added

  xix1 <- seq(1,length.out=m1$npar,by=1)
  xix2 <- seq(1+m1$npar,length.out=m2$npar,by=1)
  xinx <- seq((m1$npar+m2$npar) +1,length.out=length(unique(DivDepSpecMatrix[1:4]))-(sum(is.na(DivDepSpecMatrix[1:4]))>0) +
                length(unique(DivDepExtMatrix[1:4]))-(sum(is.na(DivDepExtMatrix[1:4]))>0),by=1)
  # So this is an array of interaction-effects as indexes into bigX. So
  # xinx <- seq(1+m1$npar+m2$npar,length.out=(4*DivDepSpecMatrix+4*DivDepExtMatrix),by=1)

  # Now let the DivDepSpecMatrix and DivDepExtMatrix be MATRICES WITH NUMBERED INDEXES. Let non-interactions be 0
  # I.e. for symmetric interactions (summed richness)
  revarix <- seq(4,length.out=sum(RE))
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

  dts = m1$dts

  myp_test <- function(x_b){
    # Or length 8 if diversity dependence in both rates internal and external
    # DIV DEP NOT IN MODEL YET:
    # I think the Ifelse won't work. Use rowSums, na.RM=T?
    lfec1 <- rowSums(cbind(drop(m1$ratefun[[1]](x_b[xix1])),
                           m1$n_norm(x_b[xix1])*x_b[xinx[DivDepSpecMatrix[1,1]]],
                           m2$n_norm(x_b[xix2])*x_b[xinx[DivDepSpecMatrix[1,2]]],
                           m1$ratefun[[6]](x_b[xix1])),na.rm=T)
    # THis means that the
    # lfec1 = base-rate from model 1 + DivDepSpecMatrix[1,1]*div1  + DivDepSpecMatrix[1,2]*div2

    lfec2 <- rowSums(cbind(drop(m2$ratefun[[1]](x_b[xix2])),
                           m1$n_norm(x_b[xix1])*x_b[xinx[DivDepSpecMatrix[2,1]]],
                           m2$n_norm(x_b[xix2])*x_b[xinx[DivDepSpecMatrix[2,2]]],
                           m2$ratefun[[6]](x_b[xix2])),na.rm=T)


    lext1 <- rowSums(cbind(drop(m1$ratefun[[2]](x_b[xix1])),
                           m1$n_norm(x_b[xix1])*x_b[xinx[DivDepExtMatrix[1,1]]],
                           m2$n_norm(x_b[xix2])*x_b[xinx[DivDepExtMatrix[1,2]]],
                           m1$ratefun[[7]](x_b[xix1])),na.rm=T)

    lext2 <- rowSums(cbind(drop(m2$ratefun[[2]](x_b[xix2])) ,
                           m1$n_norm(x_b[xix1])*x_b[xinx[DivDepExtMatrix[2,1]]],
                           m2$n_norm(x_b[xix2])*x_b[xinx[DivDepExtMatrix[2,2]]],
                           m2$ratefun[[7]](x_b[xix2])),na.rm=T)
    lsmp1 <- m1$ratefun[[3]](x_b[xix1])
    lsmp2 <- m2$ratefun[[3]](x_b[xix2])

    ll1 <- pradel_unvd_gam(ext = sapply(1:(length(dts) - 1), function(ii){
      (exp(lext1[ii]) * ((exp((exp(lfec1[ii]) - exp(lext1[ii])) * dts[ii])) - 1))/(exp(lfec1[ii]) * ((exp((exp(lfec1[ii]) - exp(lext1[ii])) * dts[ii]))) - exp(lext1[ii]))}),
      gam <- sapply(1:(length(dts) - 1), function(ii) {
        (1 - (exp(lext1[ii]) * ((exp((exp(lfec1[ii]) - exp(lext1[ii])) * dts[ii])) - 1))/(exp(lfec1[ii]) * ((exp((exp(lfec1[ii]) - exp(lext1[ii])) * dts[ii]))) - exp(lext1[ii])))/exp((exp(lfec1[ii]) - exp(lext1[ii])) * dts[ii])}),
      p = sapply(1:(length(dts)), function(ii) {rate2prob(exp(lsmp1[ii]), dts[ii])}),
      u1, n1, v1, d1)


    ll2 <- pradel_unvd_gam(ext = sapply(1:(length(dts) - 1), function(ii){
      (exp(lext2[ii]) * ((exp((exp(lfec2[ii]) - exp(lext2[ii])) * dts[ii])) - 1))/(exp(lfec2[ii]) * ((exp((exp(lfec2[ii]) - exp(lext2[ii])) * dts[ii]))) - exp(lext2[ii]))}),
      gam <- sapply(1:(length(dts) - 1), function(ii) {
        (1 - (exp(lext2[ii]) * ((exp((exp(lfec2[ii]) - exp(lext2[ii])) * dts[ii])) - 1))/(exp(lfec2[ii]) * ((exp((exp(lfec2[ii]) - exp(lext2[ii])) * dts[ii]))) - exp(lext2[ii])))/exp((exp(lfec2[ii]) - exp(lext2[ii])) * dts[ii])}),
      p = sapply(1:(length(dts)), function(ii) {rate2prob(exp(lsmp2[ii]), dts[ii])}),
      u2, n2, v2, d2)

    # Priors and hierarchies. Hplq includes both prob(RE,var) AND prior for var.
    # fph this is not correct, which(RE) gives wrong index if RE=C(F,F,T)
    hpl1 <- sapply(1:3,function(ii){ifelse(RE[ii],sum(dnorm(x_b[xix1[m1$reix[[ii]]]],0,exp(x_b[xix1[revarix[ii]]]),log=T)),0)}) #RE's m1
    hpl2 <- sapply(1:3,function(ii){ifelse(RE[ii],sum(dnorm(x_b[xix2[m2$reix[[ii]]]],0,exp(x_b[xix2[revarix[ii]]]),log=T)),0)})
    # hpl2 <- sapply(1:3,function(ii){sum(dnorm(x_b[xix2[m2$reix[[ii]]]],0,exp(x_b[xix2[3+ii]]),log=T))}) #RE's m2
    mpl1 <- sum(dnorm(x_b[xix1[1:3]],0,10),log=T) + # prior for µs
      sum(dnorm(x_b[xix2[1:3]],0,10,log=T)) + # prior for µs
      sum(dnorm(x_b[xinx],0,10,log=T)) + #prior for the interactions
      ifelse(any(RE),sum(dunif(x_b[xix1[revarix]],min=-3,max=1,log=T))+
               sum(dunif(x_b[xix2[revarix]],min=-3,max=1,log=T)),0)
    # ALSO NEED PRIORS FOR THE DRIVERS.

    p_tot = ll1$LogL + ll2$LogL + ifelse(length(hpl1)>0,sum(hpl1),0) +
      ifelse(length(hpl2)>0,sum(hpl2),0) + sum(mpl1)
    return(list(p_tot = p_tot,lraTS_1 = list(lfec=lfec1, lext=lext1,lsmp=lsmp1),
                lraTS_2 = list(lfec=lfec2,lext=lext2,lsmp=lsmp2),mpl1=mpl1,hpl1=hpl1,hpl2=hpl2,ll1=ll1,ll2=ll2))
  }
  out <- list(Clade1Mod=m1,Clade2Mod=m2,probfun=function(x){myp_test(x)$p_tot},myp_test = myp_test,npar = m1$npar+m2$npar+length(xinx),
              clade1inx = xix1, clade2inx=xix2,intinx=xinx,DivDepSpecMatrix = DivDepSpecMatrix,DivDepExtMatrix = DivDepExtMatrix)
  attr(out,"class")<-"CMR_model"
  return(out)
}


