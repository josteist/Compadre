
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
#' @param intSpec matrix of interactions between clades on speciation rates. Given as unique integers, starting at 1. For only diversity dependence of speciation rates within clades = matrix(c(1,NA,NA,2),nrow=2). For only between clade interactions = matrix(c(NA,1,2,NA),nrow=2). Similarly for intExt, starting at the ma(intSpec)+1
#' @param SpecTS1/2 drivers for speciation rates for clades 1/2
#' @param ExtTS1/2 drivers for extinction rates for clades 1/2
#' @param SmpTS1/2 drivers for sampling rates
#'
#' @return this function yields whatnow?
#' @export
#' @examples M1 <- make_BayesCMR_2clades(Obs1,Obs2,dts=rep(2,5))
#'
#'
#'
make_BayesCMR_2clades <- function(Obs1,Obs2,dts=rep(1,dim(Obs1)[2]),intSpec=matrix(rep(NA,4),nrow=2),intExt=matrix(rep(NA,4),nrow=2),
                                  SpecTS1=NULL,ExtTS1=NULL,SmpTS1=NULL,
                                  SpecTS2=NULL,ExtTS2=NULL,SmpTS2=NULL,...){
  # THis function generates a model with two interacting clades.
  # Drivers can also be fed into clades separately (SpecTS1 v SpecTS2 etc).
  # intSpec and intExt are matrices which define the interactions between clades. If NA, then no interaction. I
  # If speciation rates of clade 1 is impacte by its own and the other clades diversity (with different effects), then
  # the first row should be [1,2], with numbers indicating unique interactions. If clade1 is equally impacted in its speciation rate
  # from the richness of both own and clade 2, then same numbers can be given
  # intSpec = [1 , NA]
  #           [NA, 1]
  # denotes if interactions are between and within clades on speciation or extinction rates
  m1 <- make_BayesCMR(  Obs1, dts = dts,RE=c(T,T,T),SpecTS = SpecTS1,ExtTS= ExtTS1,SmpTS = SmpTS1,...)
  m2 <- make_BayesCMR(  Obs2, dts = dts,RE=c(T,T,T),SpecTS = SpecTS2,ExtTS= ExtTS2,SmpTS = SmpTS2,...)

  # This is the probfun.
  # ltmpfun are m1$rat
  #
  # TO DO (short term)
  # - make it allow for drivers.
  # - expand diversity effect switches to be matrices for Spec and Ext
  # How to best to that

  # - output some simple


  # TO DO (long term)
  # - allow for non-overlapping temporal spans; i.e. that one clade has 0 obs and no pars in some intervals
  # So lets just set the X's to be after eachother, first ALL the
  # x's for m1, then m2, and lastly the 1-4 pars for their interactions.
  # This will ONLY work for totally overlapping intervals of observation.
  # NOTE THAT THE MAKE_BAYESCALL was W/O diversity dependence internally, but added

  xix1 <- seq(1,length.out=m1$npar,by=1)
  xix2 <- seq(1+m1$npar,length.out=m2$npar,by=1)
  xinx <- seq((m1$npar+m2$npar) +1,length.out=length(unique(intSpec[1:4]))-(sum(is.na(intSpec[1:4]))>0) +
    length(unique(intExt[1:4]))-(sum(is.na(intExt[1:4]))>0),by=1)
  # So this is an array of interaction-effects as indexes into bigX. So
  # xinx <- seq(1+m1$npar+m2$npar,length.out=(4*intSpec+4*intExt),by=1)

 # Now let the intSpec and intExt be MATRICES WITH NUMBERED INDEXES. Let non-interactions be 0
  # I.e. for symmetric interactions (summed richness)
  myp_test <- function(x_b){
    # Or length 8 if diversity dependence in both rates internal and external
    # DIV DEP NOT IN MODEL YET:
    lfec1 <- m1$ratefunc[[1]](x_b[xix1]) + ifelse(is.na(intSpec[1,1]),0,m1$n_norm(x_b[xix1])*x_b[xinx[intSpec[1,1]]])
                                         + ifelse(is.na(intSpec[1,2]),0,m2$n_norm(x_b[xix2])*x_b[xinx[intSpec[1,2]]])
    # lfec1 = base-rate from model 1 + intSpec[1,1]*div1  + intSpec[1,2]*div2

    lfec2 <- m2$ratefunc[[1]](x_b[xix2]) + ifelse(is.na(intSpec[2,1]),0,m1$n_norm(x_b[xix1])*x_b[xinx[intSpec[2,1]]])
                                         + ifelse(is.na(intSpec[2,2]),0,m2$n_norm(x_b[xix2])*x_b[xinx[intSpec[2,2]]])


    lext1 <- m1$ratefunc[[2]](x_b[xix1]) + ifelse(is.na(intExt[1,1]),0,m1$n_norm(x_b[xix1])*x_b[xinx[intExt[1,1]]])
                                         + ifelse(is.na(intExt[1,2]),0,m2$n_norm(x_b[xix2])*x_b[xinx[intExt[1,2]]])


    lext2 <- m2$ratefunc[[2]](x_b[xix2]) + ifelse(is.na(intExt[2,1]),0,m1$n_norm(x_b[xix1])*x_b[xinx[intExt[2,1]]])
                                         + ifelse(is.na(intExt[2,2]),0,m2$n_norm(x_b[xix2])*x_b[xinx[intExt[2,2]]])
    lsmp1 <- m1$ratefunc[[3]](x_b[xix1])
    lsmp2 <- m2$ratefunc[[3]](x_b[xix2])

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


    # Priors and hierarchies.
    hpl1 <- sapply(1:3,function(ii){sum(dnorm(x_b[xix1[m1$reix[[ii]]]],0,exp(x_b[xix1[3+ii]]),log=T))}) #RE's m1
    hpl2 <- sapply(1:3,function(ii){sum(dnorm(x_b[xix2[m2$reix[[ii]]]],0,exp(x_b[xix2[3+ii]]),log=T))}) #RE's m2
    mpl1 <- sum(dnorm(x_b[xix1[1:3]],0,10),log=T) + # prior for µs
      sum(dnorm(x_b[xix2[1:3]],0,10,log=T)) + # prior for µs
      sum(dnorm(x_b[xinx],0,10,log=T)) +
      sum(dunif(x_b[xix1[4:6]],min=-3,max=1,log=T))+
      sum(dunif(x_b[xix2[4:6]],min=-3,max=1,log=T))

    p_tot = ll1$LogL + ll2$LogL + sum(hpl1) + sum(hpl2) + sum(mpl1)
    return(p_tot)
  }
  out <- list(Clade1Mod=m1,Clade2Mod=m2,probfun=myp_test,npar = m1$npar+m2$npar+length(xinx),
              clade1inx = xix1, clade2inx=xix2,intinx=xinx,intSpec = intSpec,intExt = intExt)
  attr(out,"class")<-"CMR_model"
  return(out)
}


