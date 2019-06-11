# Testing the mean scaled BayesFactor approach
# Could we extract log(BML) from this w/o model comparisons?
# FUck, er det sånn at myf_xxx functions endres når jeg redefinerer unvd?
# Yep, det ser slik ut.

 out1 <- metrop(myf_sim_fec,initial=fit$par,nbatch=10000,nspac=50,scale=stp)
 out2 <- metrop(myf_fec_ts,initial=x,nbatch=10000,nspac=50,scale=stp)

par(mfrow=c(2,1))
matplot(out1$batch,type="l")
matplot(out2$batch,type="l")
 
chain1 = out1$batch[5001:10000,];
chain2 = out2$batch[5001:10000,];

j1 <- doBayesModelLike(myf=myf_sim_fec,chain1,ndraws=1e4)
j2 <- doBayesModelLike(myf=myf_fec_ts,chain2,ndraws=1e4)

modprob1 = j1$model_prob;
normprob1 = j1$normal_prob;
use <- !(is.infinite(modprob1) | is.infinite(normprob1))
         # | is.infinite(modprob2) | is.infinite(normprob2))
mod_hat <- myf_sim_fec(j1$thetas)
norm_hat <- dmvnorm(j1$thetas,mean=j1$thetas,sigma=j1$vcv,log=T)

plot(1/seq(1,10000,by=1)*cumsum(exp(modprob1-mod_hat)/exp(normprob1-norm_hat)),type="l",xlim=c(5000,10000))
# Almost 1?

# Scale by mean(modprob1) OR modprob(mean(draws?))?
myf_sim_fec(colMeans(j1$draws))

# log(exp(normprob)/exp(modprob)) is mean(normprob1)-mean(modprob1)
logBML1 <- log(1/length(modprob1)) + (-norm_hat+mod_hat) +
  log(sum(exp(modprob1-mod_hat)/exp(normprob1-norm_hat)))



modprob2 = j2$model_prob;
normprob2 = j2$normal_prob;
use <- !(is.infinite(modprob2) | is.infinite(normprob2))
# | is.infinite(modprob2) | is.infinite(normprob2))
mod_hat2 <- myf_fec_ts(j2$thetas)
norm_hat2 <- dmvnorm(j2$thetas,mean=j2$thetas,sigma=j2$vcv,log=T)


logBML2 <- log(1/length(modprob2)) + (-norm_hat2+mod_hat2) +
  log(sum(exp(modprob2-mod_hat2)/exp(normprob2-norm_hat2)))

j <- doBayesModelLike(myf=myf_sim_fec,chain1,ndraws=100)
# Why tf doesn't this work now+
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


# >So this indicates that the TS model is has a lower BML.

# log ( sum (  ( MP(thetadraw) / MP(theta_hat) ) / (q(thetadraw)/q(theta_hat))))
# OR
# log ( sum (  ( MP(thetadraw) / MP_hat(thetadraw) ) / (q(thetadraw)/q_hat(thetadraw))))
# Think it is better with MP(theta_hat), since mean of log(probs) is not 
# good.

BF <- sum(exp(modprob1[use]-normprob1[use]-mean(modprob1[use])))/
  sum(exp(modprob2[use]-normprob2[use]-mean(modprob2[use]))) * 
  (exp(-mean(modprob2[use])+mean(modprob1[use])))



make.BayesFactor<- function(modprob1,normprob1,modprob2,normprob2){
  # Input is log model probability and log normal probability
  # with draws as in importance sampling. 
  # BF is in favour of model 1
  use <- !(is.infinite(modprob1) | is.infinite(normprob1) | is.infinite(modprob2) | is.infinite(normprob2))
  
  BF <- sum(exp(modprob1[use]-normprob1[use]-mean(modprob1[use])))/
    sum(exp(modprob2[use]-normprob2[use]-mean(modprob2[use]))) * 
    (exp(-mean(modprob2[use])+mean(modprob1[use])))
  # 
  return(BayesFactor = BF)  
}
# Could test this on a model we know better.
fakey = rnorm(25,1.2,0.3)
mod1 <- function(x){sum(dnorm(fakey,x[1],exp(x[2]),log=T))}
mod2 <- function(x){sum(dlnorm(fakey,x[1],exp(x[2]),log=T))}


x = c(0,0.1);
f1 <- optim(x,mod1,control=list(fnscale=-1))
f2 <- optim(x,mod2,control=list(fnscale=-1))
f1$value
f2$value

accs <- sapply(1:20,function(ii){metrop(mod1,initial=x,nbatch=1000,scale=0.05*ii)$accept})
stp <- which(accs<0.33)[1]*0.05;
out1 <- metrop(mod1,initial=x,nbatch=10000,nspac=50,scale=stp)

accs <- sapply(1:20,function(ii){metrop(mod2,initial=x,nbatch=1000,scale=0.05*ii)$accept})
stp <- which(accs<0.33)[1]*0.05;
out2 <- metrop(mod2,initial=x,nbatch=10000,nspac=50,scale=stp)



chain_2 =  out2$batch[round(out2$nbatch/2):out2$nbatch,];
chain_1 =  out1$batch[round(out1$nbatch/2):out1$nbatch,];
ndr = 1e6;
mbd_1 <- doBayesModelLike(mod1,chain_1,ndraws=ndr)
mbd_2 <- doBayesModelLike(mod2,chain_2,ndraws=ndr)


plot( cumsum(exp(mbd_1$model_prob-mbd_1$normal_prob))/seq(1,length(mbd_1$model_prob),by=1),type="l",ylab='')
lines(cumsum(exp(mbd_2$model_prob-mbd_2$normal_prob))/seq(1,length(mbd_2$model_prob),by=1),col='red')

bml1 <- sum(exp(mbd_1$model_prob-mbd_1$normal_prob))/ndr
bml2 <- sum(exp(mbd_2$model_prob-mbd_2$normal_prob))/ndr
bml1/bml2

bml1 <- sum(exp(mbd_1$model_prob)/exp(mbd_1$normal_prob))/ndr
bml2 <- sum(exp(mbd_2$model_prob)/exp(mbd_2$normal_prob))/ndr
bml2/bml1



modprob2  <- mbd_2$model_prob
normprob2 <- mbd_2$normal_prob
modprob1  <- mbd_1$model_prob
normprob1 <- mbd_1$normal_prob
# OK, this seems to work.
make.BayesFactor(modprob2,normprob2,modprob1,normprob1)
# Yup, they are the same.