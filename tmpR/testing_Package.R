# Testing two interacting clades.
library(Compadre)
stages = GSA_timescale[GSA_timescale$scale_level==5,]
do = seq(52,22,by=-1)
# 90 to 8 is 'Phanerozoic'
stages[do,]$interval_name
drivers <- Proxies[do[-length(do)],1:2] # not driver for last bin
dts = rev(stages$max_ma-stages$min_ma)[do]
Obs = (1*(InvertPBDB>0))[,do]
Obs <- Obs[which(rowSums(Obs)>0),]
dim(Obs)
plot(cumsum(dts),colSums(Obs>0),type="o")
mx1 <- make_BayesCMR(Obs,dts,RE=c(T,T,T),DivDep = c(F,F))
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
