# Summaries of dev feb/march 2019
# v4/5 are best, difference is how they treat the dime. v4 uses duration of
# bins to scale rates into seniorities and growth, whereas v5 uses difference
# between mid-points. They do seem to be relatively similar, but v5 perhaps
# with slightly lower bias and higher precision.


# TO make a linked model, i.e. where the diversity of one clade can be the driver of another,
# might need to augment the output of some of the model-generated stuffs. ¨¨

# For linking two clades they all must have variable rates in all three rates.
# So this options is futile. THey should be able to take same or different drivers.
# And the diversity of one is potentially a 'driver' of another. We
# get functions to get the rates $ratefunc, and could also easily make them into
# probabilities for the logLik using prade_unvd_gam.R
# So the CMR_Model should also output a function that only takes the
# rates as input (not the parameters) so that we can feed it
# rates and get a likelihood. We also need the 'full' prior to be outputted
# so we can get them independently. We already get the n_est(x) and n_est_norm(x) for
# each clade so that should be useable.

library(devtools)
setwd("\Documents")
# install.packages("Compadre")

setwd("C:/Users/josteist/Documents/Compadre")

document()
setwd("C:/Users/josteist/Documents/")
install('Compadre')


library(Compadre)
# Made a couple of changes, commit them now, push later
stages = GSA_timescale[GSA_timescale$scale_level==5,]
do = seq(90,8,by=-1)
stages[do,]$interval_name
dts = rev(stages$max_ma-stages$min_ma)[do]
Obs = (1*(InvertPBDB>0))[,do]
Obs <- Obs[which(rowSums(Obs)>0),]

m1 <- make_BayesCMR(Obs,dts=dts,RE=c(T,T,T))
m2 <- make_BayesCMR(Obs,dts=dts)
ft <- MCMC_CMR(m2,niter=1e5,nthin=200)

myESS(ft)

f1 <- MCMC_CMR(m1,niter=5e7,draweps=1e3,nthin=1e4)
f2 <- MCMC_CMR(m1,niter=1e6,draweps=1e3,nthin=1e3,cvstp = f1$Covs,
               x0=f1$Chain[5000,])
fx <- MCMC_CMR(m1,niter=1e7,draweps=1e4,nthin=1e2,vmin=1e-2)
plot(f1,botcols=(stages[do,]$color))
matplot(f1$Chain[,1:3],type="l")


library('Compadre')
S2 <- sim_BD_func_v2(spec=function(t,n){.23+sin(t)*0.1},#max(1e-8,0.8-0.1*log(n))},
                  ext = function(t,n){.2+sin(t-pi/2)*0.07},#0.4-.3*(t>12)},
                  samp = function(t,n){3.3},n_init=100,
                  dt_ints=rep(0.5,3))#rep(c(.5,4),6))
tmp <- Foote_percap(1*(S2$FosRec>0),S2$dts)
plot(cumsum(S2$dts),tmp$p_hat,type="o",ylim=c(0,.6))
lines(cumsum(S2$dts),tmp$q_hat,type="o",col='red')
lines(seq(0,sum(S2$dts),by=0.5),S2$Spec(seq(0,sum(S2$dts),by=0.5),0),lty=3)
lines(seq(0,sum(S2$dts),by=0.5),S2$Ext(seq(0,sum(S2$dts),by=0.5),0),lty=3,col='red')



S1 <- sim_BD_func_v2(spec=function(t,n){.53+sin(t)*0.3},#max(1e-8,0.8-0.1*log(n))},
                     ext = function(t,n){.4+sin(t-pi/2)*0.07},#0.4-.3*(t>12)},
                     samp = function(t,n){3.3},n_init=100,
                     dt_ints=rep(0.5,3))#rep(c(.5,4),6))
m1 <- make_BayesCMR(S1$FosRec>0,dts=S1$dts)

m2 <- make_BayesCMR_2clades(Obs1 = 1*(S1$FosRec>0),
                            Obs2 = 1*(S2$FosRec>0),dts=S2$dts,pfix=1)


f1 <- MCMC_CMR(m2,x0=runif(m2$npar,min=-.1,max=.1),niter=1e6)
matplot(f1$Chain,type="l")
matplot(f1$Chain[m2$clade1inx,],type="l")
matplot(f1$Chain[m2$clade2inx,],type="l")


# Have implemented makeCMC2clades into the package, but
# must be tested and changed so the default plot is something
# else. Could also include initial pars for sampling inside
# the makeBayes function (simpler optimization of fixed rate model)
#

m1 <- make_BayesCMR(1*(S2$FosRec),S2$dts,RE=c(T,T,T))
f1 <- MCMC_CMR(m1,niter=2e5)
m1_x <- m1;
m1_x$probfun <- function(x){m1$likfun(x)$`LogL`}
f1_x <- MCMC_CMR(m1_x,niter=2e5)
matplot(f1_x$Chain[,1:2],type="l")
matplot(f1$Chain[,1:2],type="l")



TruT <- t(sapply(1:dim(S2$Taxa)[1],
                 function(ii){S2$Taxa[ii,1]<c(cumsum(S2$dts)) & S2$Taxa[ii,2]>c(0,cumsum(S2$dts[-length(S2$dts)]))*1}))
plot(colSums(TruT),type="o",ylim=c(0,max(colSums(TruT)*1.3)))
lines(colSums(S2$FosRec>0),type="o",col='red')

# SO which one is best? For msot applications, it seems like the
# growth as immplemented in v1 is usually fine, but might give fckd results
# Change the growth from (1+l-m)^d to exp(-(l-m)*d)
# S2 <- sim_BD_dts(0.2,0.1,.2,dts = rep(c(.5,1,3,1),6))
m1 <-    make_BayesCMR(S2$FosRec>0,dts=S2$dts)
m2 <-  make_BayesCMRv2(S2$FosRec>0,dts=S2$dts)
f1 <- MCMC_CMR(m1,niter=1e4)
f2 <- MCMC_CMR(m2,niter=1e4)
# On average f4 seems worse, but the target here is moving
par(mfrow=c(1,2))
plot(f1)
plot(f2)

matplot(exp(f1$Chain[,1:3]),type="l")
matplot(exp(f2$Chain[,1:3]),type="l")
c(mean(S2$Spec(seq(0,sum(S2$dts),by=0.01),0)),
  mean(S2$Ext(seq(0,sum(S2$dts),by=0.01),0)),
  mean(S2$Samp(seq(0,sum(S2$dts),by=0.01))))


colMeans(exp(f1$Chain[-c(1:dim(f1$Chain)[1]/2),]))
colMeans(exp(f2$Chain[-c(1:dim(f1$Chain)[1]/2),]))



m1 <-    make_BayesCMR(S2$FosRec>0,dts=S2$dts,RE=c(T,T,T),DivDep=c(T,T))
m2 <-  make_BayesCMRv2(S2$FosRec>0,dts=S2$dts,RE=c(T,T,T),DivDep=c(T,T))
# It's interesting that when rates alternate clearly in sync with the intervals
# it is not detected?
x1 <-(optim(c(-1.1,-1,-1),fn=function(x){m1$likfun(c(x[1:3],rep(0,m3$npar-3)))$`LogL`},control=list(fnscale=-1))$par)
x2 <-(optim(c(-1.1,-1,-1),fn=function(x){m2$likfun(c(x[1:3],rep(0,m2$npar-3)))$`LogL`},control=list(fnscale=-1))$par)

f1a <- MCMC_CMR(m1,niter=1e6,draweps=1e3,x0=c(x1,rep(0,m1$npar-3)))
f2a <- MCMC_CMR(m2,niter=1e6,draweps=1e3,x0=c(x3,rep(0,m2$npar-3)))

  f1a <- contMCMC_CMR(f1a)
f2a <- contMCMC_CMR(f2a)

# IS it a problem if it's ill-defined if rates are exactly equal?
# well, perhaps not. It seems that 3 and 4 are closest in 'means'
tmp1 <- getRates(f1a)
tmp2 <- getRates(f2a)
c(mean(S2$Spec(seq(0,sum(S2$dts),by=0.01),0)),
  mean(S2$Ext(seq(0,sum(S2$dts),by=0.01),0)),
  mean(S2$Samp(seq(0,sum(S2$dts),by=0.01))))

par(mfrow=c(1,3))
plot(cumsum(S2$dts)[-length(S2$dts)],apply(tmp1$SpecRates,1,median),ylab='Speciation',type="o",ylim=c(0,2))
lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp2$SpecRates,1,median),type="o",col='red')
# lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp4$SpecRates,1,median),type="o",col='green')
# lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp5$SpecRates,1,median),type="o",col='yellow')
lines(seq(0,sum(S2$dts),by=0.2),S2$Spec(seq(0,sum(S2$dts),by=0.2),0),type="l",col='blue')

plot(cumsum(S2$dts),apply(tmp1$SampRates,1,median),type="o",ylab='Samp',ylim=c(0.2,.7))
lines(cumsum(S2$dts),apply(tmp2$SampRates,1,median),type="o",col='red')
abline(h=S2$Samp(1),col='blue')

plot(cumsum(S2$dts)[-length(S2$dts)],apply(tmp1$ExtRates,1,median),type="o",ylab='Ext',ylim=c(0,2))
lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp2$ExtRates,1,median),type="o",col='red')
lines(seq(0,sum(S2$dts)),S2$Ext(seq(0,sum(S2$dts)),0),type="l",col='blue')





TruT <- t(sapply(1:dim(S2$Taxa)[1],
               function(ii){S2$Taxa[ii,1]<c(cumsum(S2$dts)) & S2$Taxa[ii,2]>c(0,cumsum(S2$dts[-length(S2$dts)]))*1}))
plot(colSums(TruT),type="o",ylim=c(0,max(colSums(TruT)*1.3)))
lines(colSums(S2$FosRec>0),type="o",lty=2)

lines(colSums(S2$FosRec>0)/rate2prob(apply(tmp4$SampRates,1,median),S2$dts),type="o",col='red')
lines(colSums(S2$FosRec>0)/rate2prob(apply(tmp5$SampRates,1,median),S2$dts),type="o",col='yellow')


lines(colSums(S2$FosRec>0)/apply(tmp4$SampRates,1,median),type="o",col='red')
lines(colSums(S2$FosRec>0)/apply(tmp5$SampRates,1,median),type="o",col='red')
# So what is the 'speciation' prob:
tmp <- make_unvd(TruT)
# plot(tmp$n)
# plot(tmp$u)
# plot(tmp$v)
# plot(-(log(1-tmp$v[-1]/tmp$n[length(tmp$n)]))/S2$dts[-length(S2$dts)])
#
# plot(tmp$u[-1]/tmp$n[-21])
# plot(tmp$u[-1]/tmp$n[-21]/S2$dts[-1])
# lines(seq(0,sum(S2$dts)),S2$Spec(seq(0,sum(S2$dts)),0))
# plot(tmp$v[-1]/tmp$n[-21]/S2$dts[-1])
# lines(seq(0,sum(S2$dts)),S2$Ext(seq(0,sum(S2$dts)),0))
#
par(mfrow=c(1,2))

plot(tmp$u[-1]/tmp$n[-length(tmp$n)],ylab='new[t+1]/n[t]')
lines(rate2prob(S2$Spec(cumsum(c(0,S2$dts[-21])/2 + S2$dts/2),0),S2$dts)[-1],type="o",col='red')

plot(tmp$v/tmp$n,ylab='deat[t]/n[t]')
lines(rate2prob(S2$Ext(cumsum(c(0,S2$dts[-21])/2 + S2$dts/2),0),S2$dts),type="o",col='red')
# Here it seems like the 'data' shows higher spec rates and lower ext rates
# than the inputted rates
# Seniority is probability of a species being in i also was there
# at i-1. So 1-tmp$u/tmp$n (1 - new)
# in v1: gam = exp(-m*t)/(1 +l - m)^t)
ssp <- S2$Spec(cumsum(c(0,S2$dts[-21])/2 + S2$dts/2),0)
sex <- S2$Ext(cumsum(c(0,S2$dts[-21])/2 + S2$dts/2),0)
gam <- exp(-sex*S2$dts)/((1 + ssp - sex)^(S2$dts))
gam2 <- exp(-sex*S2$dts)/(exp(-(-ssp+sex)*S2$dts))
plot(1-tmp$u/tmp$n)
lines(gam,type="o",col='red')
  lines(gam2,type="o",col='green')
# Virtually the same SO the underestimation here is due to ones both being born and dying within and interval?
# This effectively inflates n
# So the 'bias' is dependent on the true number of species AND thus also
# the duration of the bin. As we already knew.
exp(-(-l+m)*d) #Growth from kendall
plot(f1)
#


# Can we adjust the 'weight' of some parts of the likelihood function
# to alleviate this bias?
tmp2 <- m1$likfun(x)

plot(tmp2$`Ll numerator`,tmp$n)
# obviously a clear relationship here. We want to 'weigh' each interval
# equally (or perhaps by dts) but not by number of data-points.
# THis might also improve the 'duration' problem, since larger bins
# have more species.
plot(tmp2$`Ll numerator`,
  tmp2$`Ll numerator`*(1/tmp$n/sum(tmp$n)))

sum(tmp2$`Ll numerator`)
(tmp$n/sum(tmp$n) * tmp2$`Ll numerator` )
  / sum(abs(tmp2$`Ll numerator`))



# Could we input a fit to the simulation script and simulate
# a bds process as estimated?

sum(f1a$Model$Obs[,1]) # number observed first bin, assume
