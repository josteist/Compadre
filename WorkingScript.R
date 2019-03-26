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
S1 <- sim_BD_func_v2(spec=function(t,n){0.3},#max(1e-8,0.8-0.1*log(n))},
                  ext = function(t,n){0.23},#0.4-.3*(t>12)},
                  samp = function(t,n){0.4},
                  dt_ints=rep(c(1,4),10))
m1 <- make_BayesCMR(S1$FosRec>0,dts=S1$dts)
m3 <- make_BayesCMRv3(S1$FosRec>0,dts=S1$dts)
m4 <- make_BayesCMR_v4(S1$FosRec>0,dts=S1$dts)
f1 <- MCMC_CMR(m1,niter=1e4)
f3 <- MCMC_CMR(m3,niter=1e4)
f4 <- MCMC_CMR(m4,niter=1e4)
par(mfrow=c(2,2))
plot(f1)
plot(f3)
plot(f4)
colMeans(exp(f1$Chain[501:1000,]))
colMeans(exp(f2$Chain[501:1000,]))
colMeans(exp(f3$Chain[501:1000,]))
colMeans(exp(f4$Chain[501:1000,]))

# What is a reasonable timeslice variability? From 2-12?
S1 <- sim_BD_func_v2(spec=function(t,n){0.1 + 0.002*t},
                  ext = function(t,n){0.052 + 0.0021*t},
                  dt_ints= rep(c(0.5,2,6),7),
                  samp = function(t){2})
# runif(15,min=.5,max=5))


dr1 <-
  approx(c(0,cumsum(S1$dts)),seq(0,1,length.out=(1+length(S1$dts))),
         xout = cumsum(S1$dts)/2 + c(0,cumsum(S1$dts)[-length(S1$dts)])/2)$y

m1 <- make_BayesCMR(S1$FosRec>0,  dts=S1$dts,RE=c(T,T,T),DivDep=c(F,F))
# SpecTS = t(dr1[-length(dr1)]),ExtTS=t(dr1[-length(dr1)])
m3 <- make_BayesCMRv3(S1$FosRec>0,dts=S1$dts,RE=c(T,T,T),DivDep=c(F,F))
m4 <- make_BayesCMR_v4(S1$FosRec>0,dts=S1$dts,RE=c(T,T,T),DivDep=c(F,F))
f1 <- MCMC_CMR(m1,niter=1e6)
f3 <- MCMC_CMR(m3,niter=1e6)
f4 <- MCMC_CMR(m4,niter=1e6)
par(mfrow=c(2,2))
matplot(f1$Chain[,1:3],type="l")
matplot(f3$Chain[,1:3],type="l")
matplot(f1$Chain[,7:8],type="l")
matplot(f3$Chain[,7:8],type="l")

plot(density(f1$Chain[seq(dim(f1$Chain)[1]/2,dim(f1$Chain)[1],by=10),7]))
lines(density(f1$Chain[seq(dim(f1$Chain)[1]/2,dim(f1$Chain)[1],by=10),8]),col='red')
plot( density(f3$Chain[seq(dim(f3$Chain)[1]/2,dim(f3$Chain)[1],by=10),7]))
lines(density(f3$Chain[seq(dim(f3$Chain)[1]/2,dim(f3$Chain)[1],by=10),8]),col='red')



m1$probfun(x)
m4$probfun(x)

plot(f1,log=F)
plot(f3)



?MCMC_CMR


setwd("C:/Users/josteist/Documents/Compadre/R")
