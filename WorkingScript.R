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


setwd("Documents")
install.packages("Compadre")
setwd("./Compadre")
document()


# Made a couple of changes, commit them now, push later



S1 <- sim_BD_func(spec=function(t,n){0.8},
                  ext = function(t,n){0.64},
                  dt_ints=rep(2,12))
m1 <- make_BayesCMR(S1$FosRec>0,dts=S1$dts)
m2 <- make_BayesCMRv2(S1$FosRec>0,dts=S1$dts)
m3 <- make_BayesCMRv3(S1$FosRec>0,dts=S1$dts)
m4 <- make_BayesCMR_old(S1$FosRec>0,dts=S1$dts)
f1 <- MCMC_CMR(m1,niter=1e4)
f2 <- MCMC_CMR(m2,niter=1e4)
f3 <- MCMC_CMR(m3,niter=1e4)
f4 <- MCMC_CMR(m4,niter=1e4)
par(mfrow=c(2,2))
plot(f1)
plot(f2)
plot(f3)
plot(f4)


m1 <- make_BayesCMR(S1$FosRec>0,dts=S1$dts,RE=c(T,T,T))
m3 <- make_BayesCMRv3(S1$FosRec>0,dts=S1$dts,RE=c(T,T,T))
f1 <- MCMC_CMR(m1,niter=2e5)
f3 <- MCMC_CMR(m3,niter=2e5)
par(mfrow=c(2,2))
plot(f1)
plot(f3)



?MCMC_CMR
