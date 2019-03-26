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
setwd("C:/Users/jostest/Documents/Compadre")
document()
setwd("C:/Users/josteist/Documents/")
install('Compadre')


library(Compadre)
# Made a couple of changes, commit them now, push later


library('Compadre')
S2 <- sim_BD_func_v2(spec=function(t,n){.5+sin(t)*0.4},#max(1e-8,0.8-0.1*log(n))},
                  ext = function(t,n){.4+cos(t-pi/2)*0.3},#0.4-.3*(t>12)},
                  samp = function(t,n){0.4},
                  dt_ints=rep(c(.5,2,3,1),5),n_init=200)
TruT <- t(sapply(1:dim(S2$Taxa)[1],
                 function(ii){S2$Taxa[ii,1]<c(cumsum(S2$dts)) & S2$Taxa[ii,2]>c(0,cumsum(S2$dts[-length(S2$dts)]))*1}))
plot(colSums(TruT),type="o",ylim=c(0,max(colSums(TruT)*1.3)))
lines(colSums(S2$FosRec>0),type="o",col='red')

# SO which one is best? For msot applications, it seems like the
# growth as immplemented in v1 is usually fine, but might give fckd results
# Change the growth from (1+l-m)^d to exp(-(l-m)*d)
# S2 <- sim_BD_dts(0.2,0.1,.2,dts = rep(c(.5,1,3,1),6))
m1 <-    make_BayesCMR(S2$FosRec>0,dts=S2$dts)
m3 <-  make_BayesCMRv3(S2$FosRec>0,dts=S2$dts)
m4 <- make_BayesCMR_v4(S2$FosRec>0,dts=S2$dts)
m5 <- make_BayesCMR_v5(S2$FosRec>0,dts=S2$dts)
f1 <- MCMC_CMR(m1,niter=3e4)
f3 <- MCMC_CMR(m3,niter=3e4)
f4 <- MCMC_CMR(m4,niter=3e4)
f5 <- MCMC_CMR(m4,niter=3e4)
# On average f4 seems worse, but the target here is moving
par(mfrow=c(2,2))
plot(f1)
plot(f3)
plot(f4)
plot(f5)

par(mfrow=c(2,2))
matplot(exp(f1$Chain[,1:3]),type="l")
matplot(exp(f3$Chain[,1:3]),type="l")
matplot(exp(f4$Chain[,1:3]),type="l")
matplot(exp(f5$Chain[,1:3]),type="l")
c(mean(S2$Spec(seq(0,sum(S2$dts),by=0.01),0)),
  mean(S2$Ext(seq(0,sum(S2$dts),by=0.01),0)),
  mean(S2$Samp(seq(0,sum(S2$dts),by=0.01))))


  par(mfrow=c(1,1))
  colMeans(exp(f1$Chain[-c(1:dim(f1$Chain)[1]/2),]))
  colMeans(exp(f3$Chain[-c(1:dim(f1$Chain)[1]/2),]))
  colMeans(exp(f4$Chain[-c(1:dim(f1$Chain)[1]/2),]))
  colMeans(exp(f5$Chain[-c(1:dim(f1$Chain)[1]/2),]))

# HEre v4 and 5 are obviously better. But all underestimate the sampling.
# What is a reasonable timeslice variability? From 2-12?
# S2 <- sim_BD_func_v2(spec=function(t,n){0.3 + 0.1*t},#((t %% 2)>1)*0.25},
#                   ext = function(t,n){0.2 + 0.1*t},
#                   dt_ints= rep(c(.3,3),10),
#                   samp = function(t){1-0.015*t})
# FOr non-varying interval durs, 1&3 seems to do a bit better, but all
  # runif(15,min=.5,max=5))
# v3 seems best, much less 'variability'.

S2 <- sim_BD_func_v2(spec=function(t,n){0.5+0.8*(t>22)},#max(1e-8,0.8-0.1*log(n))},
                     ext = function(t,n){.45+0.95*(t>22)},#0.4-.3*(t>12)},
                     samp = function(t,n){0.4},
                     dt_ints=rep(c(.5,2,3,1),5),n_init=200)


m1 <-    make_BayesCMR(S2$FosRec>0,dts=S2$dts,RE=c(T,T,T))
m3 <-  make_BayesCMRv3(S2$FosRec>0,dts=S2$dts,RE=c(T,T,T))
m4 <- make_BayesCMR_v4(S2$FosRec>0,dts=S2$dts,RE=c(T,T,T))
m5 <- make_BayesCMR_v5(S2$FosRec>0,dts=S2$dts,RE=c(T,T,T))
# It's interesting that when rates alternate clearly in sync with the intervals
# it is not detected?
x1 <-(optim(c(-1.1,-1,-1),fn=function(x){m1$likfun(c(x[1:3],rep(0,m3$npar-3)))$`LogL`},control=list(fnscale=-1))$par)
x3 <-(optim(c(-1.1,-1,-1),fn=function(x){m3$likfun(c(x[1:3],rep(0,m3$npar-3)))$`LogL`},control=list(fnscale=-1))$par)
x4 <-(optim(c(-1.1,-1,-1),fn=function(x){m4$likfun(c(x[1:3],rep(0,m3$npar-3)))$`LogL`},control=list(fnscale=-1))$par)
x5 <-(optim(c(-1.1,-1,-1),fn=function(x){m5$likfun(c(x[1:3],rep(0,m3$npar-3)))$`LogL`},control=list(fnscale=-1))$par)

f1a <- MCMC_CMR(m1,niter=5e5,draweps=1e4,x0=c(x1,rep(0,m1$npar-3)))
f3a <- MCMC_CMR(m3,niter=5e5,draweps=1e4,x0=c(x3,rep(0,m1$npar-3)))
f4a <- MCMC_CMR(m4,niter=5e5,draweps=1e4,x0=c(x4,rep(0,m1$npar-3)))
f5a <- MCMC_CMR(m5,niter=5e5,draweps=1e4,x0=c(x5,rep(0,m1$npar-3)))

# IS it a problem if it's ill-defined if rates are exactly equal?
# well, perhaps not. It seems that 3 and 4 are closest in 'means'
f1a <- MCMC_CMR(m1,niter=5e5,cvstp = f1a$Covs,draweps=5e5/100,x0=c(x1,rep(0,m1$npar-3)))
f3a <- MCMC_CMR(m3,niter=5e5,cvstp = f3a$Covs,draweps=5e5/100,x0=c(x3,rep(0,m1$npar-3)))
f4a <- MCMC_CMR(m4,niter=5e5,cvstp = f4a$Covs,draweps=5e5/100,x0=c(x4,rep(0,m1$npar-3)))
f5a <- MCMC_CMR(m5,niter=5e5,cvstp = f5a$Covs,draweps=5e5/100,x0=c(x5,rep(0,m1$npar-3)))

tmp1 <- getRates(f1a)
tmp3 <- getRates(f3a)
tmp4 <- getRates(f4a)
tmp5 <- getRates(f5a)
c(mean(S2$Spec(seq(0,sum(S2$dts),by=0.01),0)),
  mean(S2$Ext(seq(0,sum(S2$dts),by=0.01),0)),
  mean(S2$Samp(seq(0,sum(S2$dts),by=0.01))))

par(mfrow=c(1,3))
plot(cumsum(S2$dts)[-length(S2$dts)],apply(tmp1$SpecRates,1,median),ylab='Speciation',type="o",ylim=c(0,2))
lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp3$SpecRates,1,median),type="o",col='red')
lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp4$SpecRates,1,median),type="o",col='green')
lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp5$SpecRates,1,median),type="o",col='yellow')
lines(seq(0,sum(S2$dts),by=0.2),S2$Spec(seq(0,sum(S2$dts),by=0.2),0),type="l",col='blue')

plot(cumsum(S2$dts),apply(tmp1$SampRates,1,median),type="o",ylab='Samp',ylim=c(0.2,.7))
lines(cumsum(S2$dts),apply(tmp3$SampRates,1,median),type="o",col='red')
lines(cumsum(S2$dts),apply(tmp4$SampRates,1,median),type="o",col='green')
lines(cumsum(S2$dts),apply(tmp5$SampRates,1,median),type="o",col='yellow')
# lines(seq(0,sum(S2$dts),by=0.1),S2$Samp(seq(0,sum(S2$dts),by=0.1)),type="l",col='blue')
abline(h=S2$Samp(1),col='blue')

plot(cumsum(S2$dts)[-length(S2$dts)],apply(tmp1$ExtRates,1,median),type="o",ylab='Ext',ylim=c(0,2))
lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp3$ExtRates,1,median),type="o",col='red')
lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp4$ExtRates,1,median),type="o",col='green')
lines(cumsum(S2$dts)[-length(S2$dts)],apply(tmp5$ExtRates,1,median),type="o",col='yellow')
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


