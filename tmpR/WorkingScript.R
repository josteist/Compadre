library(devtools)
setwd('C:/Users/josteist/Documents/Compadre')
document()
setwd('C:/Users/josteist/Documents')
install('Compadre')

library(Compadre)

# Testing two interacting clades.
rm(list=ls())
library(Compadre)
set.seed(211080+2)
stages = GSA_timescale[GSA_timescale$scale_level==5,]
do = seq(90,78,by=-1)
# 90 to 8 is 'Phanerozoic'
stages[do,]$interval_name
drivers <- Proxies[do,] # not driver for last bin
dts = stages[do,]$max_ma-stages[do,]$min_ma
names(dts) <- stages[do,]$interval_name


Obs = (t(Occ_species>0))[,do]
Obs <- Obs[which(rowSums(Obs)>0),]

dim(Obs)
m1 <- make_BayesCMR(Obs,dts,spec=~div*SeaLev,ext = ~time,samp = ~time,data=drivers)

x = runif(m1$npar)
m1$mmspec(x) %*% x[m1$inx$specInx]
m1$specfun(x)
# This seems to work yes...

f1 <- MCMC_CMR(m1,niter=5e5,nthin=1,draweps=1e4,vmin=1e-3)
# matplot(f1$Chain[,1:3],type="l")
# matplot(f1$Chain[,31:33],type="l")

vmin = 1e-2
(2.38/(sqrt(m1$npar)))^2*vmin
hist(diag(f1$Covs),xlim=c(0,vmin))
points(0,(2.38/(sqrt(m1$npar)))^2*vmin,col='red',pch=20)
(2.38/(sqrt(m1$npar)))^2*vmin
# save.image('LongRun030819.RData')
# Even this did not converge properly. It took ~3.5 days.


# Testing the mcmc package.
install.packages('mcmc')
library(mcmc)
vmin = 0.05
f <- metrop(m1$probfun,initial = c(-2,-2.1,-1.9,runif(min=-0.01,max=0.01,m1$npar-3)),
            nbatch=1e5,nspac=10,scale=(2.38/(sqrt(m1$npar)))^2*vmin)
f$accept
vmins <- seq(1e-4,1e-1,length.out=10);

tmp <- sapply(1:length(vmins),function(ii){
  metrop(m1$probfun,initial = c(-2,-2.1,-1.9,runif(min=-0.01,max=0.01,m1$npar-3)),
            nbatch=1e3,nspac=1,scale=(2.38/(sqrt(m1$npar)))^2*vmins[ii])$accept})
plot(vmins,tmp)
vmins[4]
vmin = 0.05
f <- metrop(m1$probfun,initial = c(-2,-2.1,-1.9,runif(min=-0.01,max=0.01,m1$npar-3)),
            nbatch=1e5,nspac=10,scale=(2.38/(sqrt(m1$npar)))^2*vmins[4])

cvstp <- (2.38/(sqrt(m1$npar)))^2*(cov(f$batch)  +
                                           diag(vmin,m1$npar))
f2 <- metrop(m1$probfun,initial = c(-2,-2.1,-1.9,runif(min=-0.01,max=0.01,m1$npar-3)),
            nbatch=1e4,nspac=1,scale=cvstp)

matplot(f2$batch[,1:3],type="l")
matplot(f2$batch[,34:36],type="l")


f1$Chain = f$batch
plot(f1)
coda::effectiveSize(f$batch)


m1 <- make_BayesCMR(Obs,dts,spec=~time+div*SeaLev,ext = ~time,samp = ~time,data=drivers)
m1$specfun
m1$inx$specInx
x = runif(m1$npar)
x[6] =0
# m1$mmspec(x)
tmp1 <- cbind(m1$mmspec(x) %*% x[m1$inx$specInx])
x[6] =1
tmp2 <- cbind(m1$mmspec(x) %*% x[m1$inx$specInx])
plot(drivers$SeaLev[-length(dts)],tmp2)
points(drivers$SeaLev[-length(dts)],tmp1,col='red')

f1 <- MCMC_CMR(m1,niter=5e5,nthin=1,draweps=1e4,vmin=1e-3)



Obs1 = (t(Occ_genera>0))[,do]
Obs1 <- Obs1[which(rowSums(Obs1)>0),]

Obs1 = (t(Occ_genera[,which(Taxonomy_genera$phylum == 'Brachiopoda')]>0))[,do]
Obs1 <- Obs1[which(rowSums(Obs1)>0),]


Obs2 = (t(Occ_species>0))[,do]
Obs2 <- Obs2[which(rowSums(Obs2)>0),]

dim(Obs)
m1 <- make_BayesCMR(Obs1,dts,spec=~time+div,ext = ~time,samp = ~time)
m2 <- make_BayesCMR(Obs2,dts)#,spec=~time,ext = ~time,samp = ~time)
f1 <- MCMC_CMR(m1,niter=1e4)
f2 <- MCMC_CMR(m2,niter=1e4)
plot(f1,stages = stages[do,])
plot(f2,stages = stages[do,])

ma1 <- make_BayesCMR(Obs,dts,  spec= ~ d180_cor+time,ext = ~ d180_cor*d13C*div+time,samp = ~time,data = drivers)
ma2 <- make_BayesCMR(Obs,dts,  spec= ~ d180_cor+time,ext = ~ d180_cor*d13C*div+time,samp = ~time,data = drivers)

fa1 <- MCMC_CMR(ma1,niter=1e6)
fa2 <- MCMC_CMR(ma2,niter=1e6)



# Doing taxonomic subsets.
Obs1 = (t(Occ_genera[,which(Taxonomy_genera$phylum == 'Brachiopoda')]>0))[,do]
Obs1 <- Obs1[which(rowSums(Obs1)>0),]

Obs2 = (t(Occ_genera[,which(Taxonomy_genera$class == 'Bivalvia')]>0))[,do]
Obs2 <- Obs2[which(rowSums(Obs2)>0),]

plot(rev(stages[do,]$max_ma/2+stages[do,]$min_ma/2),
     colSums(Obs1>0),type="o",xlim=rev(range(stages[do,]$max_ma)))
lines(rev(stages[do,]$max_ma/2+stages[do,]$min_ma/2),
     colSums(Obs2>0),type="o",xlim=rev(range(stages[do,]$max_ma)),col='red')


m1 <- make_BayesCMR(Obs1,dts,spec=~time, ext=~time,samp=~time)
m2 <- make_BayesCMR(Obs2,dts,spec=~time, ext=~time,samp=~time)


f1 <- MCMC_CMR(m1,niter=1e6,nthin=1e3)
f2 <- MCMC_CMR(m2,niter=1e6,nthin=1e3)

fax <- contMCMC_CMR(fa)

plotDrivers(fa,samp=F,nsamp=500)

par(mfrow=c(2,4))
for (ii in 1:length(ma$inx$extInx)){
hist(fa$Chain[-c(1:dim(fa$Chain)[1]/2),ma$inx$extInx[ii]],main=names(ma$inx$extInx)[ii])
abline(v=0,col='red')
}

par(mfrow=c(2,3))
for (ii in 1:length(ma$inx$specInx)){
  hist(fa$Chain[-c(1:dim(fa$Chain)[1]/2),ma$inx$specInx[ii]],main=names(ma$inx$specInx)[ii])
  abline(v=0,col='red')
}

smp <- seq(dim(fa$Chain)[1]/2,dim(fa$Chain)[1],by=100)



rowSums(model.matrix(~d180_cor+div*d13C,drivers) %*% x[ma$inx$specInx])
ma$specfun(x)
