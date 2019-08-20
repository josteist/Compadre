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
make_BayesCMR_2cladd <- function(Obs1,Obs2,dts=rep(1,dim(Obs1)[2]),
                                 spec1 = ~time, ext1 = ~time, samp1 = ~time,
                                 spec2 = ~time, ext2 = ~time, samp2 = ~time,data=NULL,...){
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


  # Perhaps call this w/o any formula for spec/ext, and implement those outside?
  # OR actually feed in data with div1 and div2 (then treated as 'drivers'), and repopulate the data.frames outside?
  data_tmp <- data.frame(cbind(drivers,div1=rep(0.1,length(dts)),div2=rep(0.2,length(dts)),div12 = rep(0.3,length(dts))))
  spec1 <- ~ div1*d180_cor + d13C
  ext1   <- ~time + div2*d13C
  samp1 <- ~time
  m1 <- make_BayesCMR(  Obs1, dts = dts,spec=spec1,ext=ext1,samp=samp1,data=data_tmp)
  m2 <- make_BayesCMR(  Obs2, dts = dts,spec=spec2,ext=ext2,samp=samp2,data=data_tmp)

  # This is the probfun.
  # ltmpfun are m1$rat
  #

  bigdrivers <- function(x){data.frame(drivers,div1=m1$n_norm(x[xinx1]),div2=m2$n_norm(x[xinx2]),
                                       div12=m1$n_norm(x[xinx1])*m2$n_norm(x[xinx2]))}


  rowSums(model.matrix(spec1,data_tmp) %*% x[ma$inx$specInx])
  ma$specfun(x)

  model.matrix(~1+d13C*div1 + div2,data=drivers)


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
  revarix <- rep(NA,3);
  revarix[RE] = seq(4,length.out=sum(RE));
  # seq(4,length.out=sum(RE))
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
      ifelse(any(RE),sum(dunif(x_b[xix1[revarix]],min=-3,max=1,log=T),na.rm=T)+
               sum(dunif(x_b[xix2[revarix]],min=-3,max=1,log=T),na.rm=T),0)
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


