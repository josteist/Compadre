# Testing two interacting clades.
library(Compadre)
stages = GSA_timescale[GSA_timescale$scale_level==5,]
do = seq(90,8,by=-1)
# 90 to 8 is 'Phanerozoic'
stages[do,]$interval_name
drivers <- Proxies[do[-length(do)],1:2] # not driver for last bin
dts = rev(stages$max_ma-stages$min_ma)[do]
Obs = (1*(InvertPBDB>0))[,do]
Obs <- Obs[which(rowSums(Obs)>0),]
dim(Obs)


plot(cumsum(dts),colSums(Obs>0),type="o")
set.seed(190480)

# Splitting randomly
O1 <- Obs[1:20000,];
O2 <- Obs[20001:61834,]
dim(O1)

M1 <- make_BayesCMR_2clades(O1,O2,dts = dts,RE=c(T,T,T),
                            DivDepSpecMatrix = matrix(c(1,2,1,2),ncol=2),
                            DivDepExtMatrix = matrix(c(3,4,3,4),ncol=2),
                            SpecTS1 = t(drivers),
                            ExtTS2  = t(drivers[,1]),
                            Driv_x_Div_Spec_1 = c(T,T),
                            Driv_x_Div_Ext_2  = c(T))
xx = runif(M1$npar,min=-0.01,max=0.1)
xx[M1$intinx]=0.1
M1$myp_test(xx)

xx = runif(M1$npar,min=-1,max=1)
M1$myp_test(xx)$ll1
# It depends on the parameters, how the fuck did that come about?
# SO I guess it's one of the rates that are pushing against the boundary.
# gam is almost 0 before it fails, and xi is almost 1.
# This might be prevented with initials around 0
M1$myp_test(xx)$lrats2

M1$probfun(xx)


M1$Clade1Mod$specfun(xx[M1$clade1inx])
M1$Clade1Mod$ratefun[[2]](xx[M1$clade1inx])

hpl1 <- sapply(1:3,function(ii){
  ifelse(RE[ii],sum(dnorm(xx[M1$clade1inx[M1$Clade1Mod$reix[[ii]]]],0,
                      exp(xx[M1$Clade1Mod$reix[[ii]]]),log=T)),0)}) #RE's m1

# intExt = matrix(c(3,4,3,4),ncol=2),#,
                            # SpecTS1 = drivers[,1],SpecTS2=drivers[,2],
                            # SpecInt1 = c(T))#,ExtInt2=c(T))
# xx <- runif(M1$npar)
# M1$Clade1Mod$ratefun[[1]](xx)
# M1$probfun(xx)
ma <- make_BayesCMR(O1,dts=dts,DivDep=c(T,T),SpecTS = t(drivers),SpecInt =c(T,T))
mb <- make_BayesCMR(O2,dts=dts,DivDep=c(T,T),ExtTS = t(drivers[,1]),ExtInt=c(T))
# Why is this Soooooo slow?
fx1 <- MCMC_CMR(M1,vmin=1e-5,x0=runif(M1$npar,min=-0.05,max=0.05),
                niter=1e5,nthin=100)

# While slightly strange, all interaction coeffs are >0, i.e.
# positive, i.e. more speciation with more species....
# Error in cbind(m2$ratefun[[2]](x_b[xix2]), m1$n_norm(x_b[xix1]) *
# x_b[xinx[intExt[2,  :
# number of rows of matrices must match (see arg 4)

fxa <- MCMC_CMR(ma,vmin=1e-4,niter=5e4)
fxb <- MCMC_CMR(mb,vmin=1e-4,niter=5e4)
par(mfrow=c(1,4))
matplot(fx1$Chain[-c(1:300),c(M1$clade1inx,M1$intinx[c(1,3)])],type="l")
matplot(fxa$Chain[,c(1:3,6,7,8,9,4,5)],type="l")
matplot(fx1$Chain[-c(1:300),c(M1$clade2inx,M1$intinx[c(2,4)])],type="l")
matplot(fxb$Chain[,c(1:3,6,7,4,5)],type="l")
# Seems to fit w/o div_driver interactions, also with driver interactions
# they are the same. Testing for interaction and interactions complete.

par(mfrow=c(1,4))
matplot(fx1$Chain[-c(1:300),c(7,9)],type="l")
matplot(fxa$Chain[-c(1:300),c(4,5)],type="l")
matplot(fx1$Chain[-c(1:300),c(8,11)],type="l")
matplot(fxb$Chain[-c(1:300),c(4,5)],type="l")

# So here, it seems like we have both rates positively diversity
# dependent? That's raaaaather strange. Anyhows,
# Whan rates are simply divdep within each clade, results are
# same with the clade2func. Now testing with all impacts on
# clade 1 (intSpec=c(1,1,NA,2)),intExt(c(3,4,NA,5)) and only
# self-impact of clade 2
cbind(M1$Clade1Mod$ratefun[[2]](xx[M1$clade1inx]),
      M1$Clade1Mod$n_norm(xx[M1$clade1inx])*xx[M1$intinx[intExt[1,1]]],
      M1$Clade2Mod$n_norm(xx[M1$clade2inx])*xx[M1$intinx[intExt[1,2]]],
      M1$Clade1Mod$ratefun[[7]](xx[M1$clade1inx]))


cbind(drop(M1$Clade2Mod$ratefun[[2]](xx[M1$clade2inx])),
      M1$Clade1Mod$n_norm(xx[M1$clade1inx])*xx[M1$intinx[intExt[2,1]]],
      M1$Clade2Mod$n_norm(xx[M1$clade2inx])*xx[M1$intinx[intExt[2,2]]],
      M1$Clade2Mod$ratefun[[7]](xx[M2$clade1inx]))


mx1 <- make_BayesCMR(Obs,dts,RE=c(T,F,T),DivDep = c(F,F))
x0  <- getInit(mx1)
ESS(fx1)
matplot(fx1$Chain,type="l")
# I NEED TO CHANGE DEFAULT vmin; leads to no change for simple 3-rate models
fx1a <- MCMC_CMR(mx1,x0=x0,vmin=1e-5)
fx1b <- MCMC_CMR(mx1,niter = 1e5,draweps=500,vmin=1e-5,
                 x0=x0+runif(length(x0),min=-.1,max=.1))

matplot(fx1$Chain[1:1e4,1:3],type="l")
matplot(fx1b$Chain[,1:3],type="l")



x <- fx1b$Chain[1,]
x2 <-   sapply(7:mx1$npar,
         function(jj){optim(c(-1),function(y){-mx1$probfun(c(x[1:jj],y,x[(jj+1):mx1$npar]))},method='Brent',lower=-10,upper=10)$par})

plot(fx1)
ESS(fx1)

mx2 <- make_BayesCMR(Obs,dts,RE=c(T,T,T),DivDep = c(T,F))
mx3 <- make_BayesCMR(Obs,dts,RE=c(T,T,T),DivDep = c(F,T))
mx4 <- make_BayesCMR(Obs,dts,RE=c(T,T,T),DivDep = c(T,T))
mx5 <- make_BayesCMR(Obs,dts,RE=c(T,T,T),DivDep = c(F,F),SpecTS = t(drivers),ExtTS=t(drivers))

set.seed(190480)



fx1  <- MCMC_CMR(mx1,niter=1e5,nthin=10,draweps=200,vmin=1e-2)
fx1a <- MCMC_CMR(mx1)

fx1b <- MCMC_CMR_tmp(mx1,niter=4e5,nthin=10,draweps=1000,vmin=1e-2)

  fx1b <- MCMC_CMR(mx1,niter=1e5,nthin=10,draweps=100,vmin=1e-3)
fx1c <- MCMC_CMR(mx1,niter=1e5,nthin=10,draweps=1000,vmin=1e-3)

(2.38/(sqrt(mx1$npar)))^2 *1e-3
hist(diag(fx1c$Covs),21)

fx1a <- contMCMC_CMR(fx1)
fx2 <- MCMC_CMR(mx2,niter=5e5,nthin=50)
fx3 <- MCMC_CMR(mx3,niter=5e5,nthin=50)
fx4 <- MCMC_CMR(mx4,niter=5e5,nthin=50)
fx5 <- MCMC_CMR(mx5,niter=5e5,nthin=50)


# Seems like MASS does not have the dmvnom
## Testing model likelihood
doBayesModelLike <- function(myf,chain,ndraws=1e4){
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

bml1 <- replicate(20,
                  doBayesModelLike(mx1$probfun,fx1$Chain[seq(5000,10000,by=5),],ndraws=5000)$logBML)
bml2 <- replicate(20,
                  doBayesModelLike(mx2$probfun,fx2$Chain[seq(5000,10000,by=5),],ndraws=5000)$logBML)
bml3 <- replicate(20,
                  doBayesModelLike(mx3$probfun,fx3$Chain[seq(5000,10000,by=5),],ndraws=5000)$logBML)
bml4 <- replicate(20,
                  doBayesModelLike(mx4$probfun,fx4$Chain[seq(5000,10000,by=5),],ndraws=5000)$logBML)
bml5 <- replicate(20,
                  doBayesModelLike(mx5$probfun,fx5$Chain[seq(5000,10000,by=5),],ndraws=5000)$logBML)

boxplot(cbind(bml1,bml2,bml3,bml4,bml5))
# Even this seems to work relatively well now. Do the drivers alone?
# For this case, the difference in bml is not too big.
Obs1 <- Obs[1:30000,30:43];
Obs2 <- Obs[30001:60000,30:43];
dtsx <- dts[30:43]
Obs1 <- Obs1[which(rowSums(Obs1)>0),]
Obs2 <- Obs2[which(rowSums(Obs2)>0),]

intS <- matrix(c(1,NA,NA,2),nrow=2)
intE <- matrix(c(3,5,4,6),nrow=2)
mB <- make_BayesCMR_2clades_tmp(Obs1,Obs2,dts=dtsx,intSpec=intS,intExt = intE)

fB <- MCMC_CMR(mB,niter=1e4)

m1 <- make_BayesCMR(Obs1,dtsx,RE=c(T,T,T))
m2 <- make_BayesCMR(Obs2,dtsx,RE=c(T,T,T))


f1 <- MCMC_CMR(m1,niter=4e5)#,x=runif(m1$npar,min=-0.01,max=0.02))
f1a <- contMCMC_CMR(f1,niter=1e4);

ESS(f2)
f2 <- MCMC_CMR(m2,niter=5e4)


mB <- make_BayesCMR_2clades(Obs1,Obs2,dts=dtsx,
                            SpecTS1 = t(runif(length(dtsx)-1,min=-1,max=2)),
                            ExtTS2  = t(runif(length(dtsx)-1,min=0.1,max=0.3)))

fB <- MCMC_CMR(mB,niter=1e3)

#
#
#
# myESS <- function(f1){
#   coda::effectiveSize(as.mcmc(f1$Chain[,-c(1:dim(f1$Chain)[2]/2)]))
# }



stages = GSA_timescale[GSA_timescale$scale_level==5,]
mp = stages$min_ma/2 + stages$max_ma/2; # midpoints of stages.

# As example we will analyse the invertbare fauna leading up to the Permian extinction.
use = seq(61,49,by=-1); # these are the stages in the whole of Permian and firt half Triassic.
Ocs <- InvertPBDB[,use]; # Taking out these stages
Ocs <- Ocs[rowSums(Ocs)>0,]; # keeping only species with obs here.
set.seed(190480)
#
plot(mp[use],colSums(Ocs>0),type="o",xlim=rev(range(mp[use])),log='',ylab='# species')
dts <- stages$max_ma[use] - stages$min_ma[use]
cls <- stages$color[use];
# Simple model
m1 <- make_BayesCMR(Ocs,dts)
f1 <- MCMC_CMR(m1,vmin=1e-4)
ESS(f1) # check effective sample size.
matplot(f1$Chain,type="l")
plot(f1)
# So, how 'large' was the mass extinction. The transition out from the Permian
# is the 9th transition here. Making a 'mock'-driver with 0's all other transitions
# than this.
driver = rep(0,length(dts)-1); # All 0. Drivers have a length of length(dts)-1, since there is not transition out of the last interval.
driver[9] = 1;
m2 <- make_BayesCMR(Ocs,dts,ExtTS = driver)
# Check is the driver is implementer
m2

f2 <- MCMC_CMR(m2,vmin=1e-4,niter=5e4)
ESS(f2)
tmp <- plot(f2,log=T,botcol=cls,max_ma=max(stages[use,]$max_ma))
# Getting median extinction rates.
apply(tmp$ExtRates,1,median)
# Calculting the extinction probability for this transition
rate2prob(median(tmp$ExtRates[9,]),dts[9])
# Actually, we should perhaps rather use the 'Raup' stuff here too, instead of simple rate2prob?

# Trying a similar model, but now including temporally varying sampling rates. Is the imapct of mass extinction exacerbated by ignoring different sampling rates over time?

# RE consists of three TRUE/FALSE statements on whether or not to
# include random effects/time-varying rates for Speciation, extinction and/or Sampling rates.
# We want to only include for sampling rates, so we use c(F,F,T))
m3 <- make_BayesCMR(Ocs,dts,ExtTS = driver,RE=c(F,F,T))
# Check is the driver is implementer
m3
f3 <- MCMC_CMR(m3,vmin=1e-4,niter=1e5)
tmp <- plot(f3,log=T,botcol=cls,max_ma=max(stages[use,]$max_ma))
# Getting median extinction rates.
apply(tmp$ExtRates,1,median)
# A much lower extinction rate. Seems like the impact of the Permian on
# the species level here might have been exaggerated due to assumed homogenous sampling rate.


# Final big model with all temporally varying drivers. This can be implemented
# in two different way; either with random effects for all transitions. Then, on
# the log scale, all deviations from the global extinction rate are assumed to be
# normally distributed with a common variance. If we implement a 'driver' with
# only 1 at the focal transition, then this transition is estimated a separate parameters
# and thus not part of the 'random' effects. In principle it's then a fixed effect with it's own prior.

m4a <- make_BayesCMR(Ocs,dts,RE=c(T,T,T))
m4b <- make_BayesCMR(Ocs,dts,RE=c(T,T,T),ExtTS = driver,ExtInt=c(T))
m4b

f4a <- MCMC_CMR(m4a,vmin=1e-4,niter=1e6)
f4b <- MCMC_CMR(m4b,vmin=1e-4,niter=1e6)
matplot(f4b$Chain[,c(1:3,7)],type="l")
#

plot(f4a,log=T,botcols=cls,max_ma=max(stages[use,]$max_ma))
plot(f4b,log=T,botcols=cls,max_ma=max(stages[use,]$max_ma))
