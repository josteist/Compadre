# Two clade tests
rm(list=ls())
library(Compadre)
do = seq(50,8,by=-1)
# 90 to 8 is 'Phanerozoic'
stages = GSA_timescale[GSA_timescale$scale_level==5,]

stages[do,]$interval_name
drivers <- Proxies[do,] # not driver for last bin
dts = stages[do,]$max_ma-stages[do,]$min_ma
names(dts) <- stages[do,]$interval_name

# Doing taxonomic subsets.
Obs1 = (t(Occ_genera[,which(Taxonomy_genera$phylum == 'Brachiopoda')]>0))[,do]
Obs1 <- Obs1[which(rowSums(Obs1)>0),]

Obs2 = (t(Occ_genera[,which(Taxonomy_genera$class == 'Bivalvia')]>0))[,do]
Obs2 <- Obs2[which(rowSums(Obs2)>0),]

plot(rev(stages[do,]$max_ma/2+stages[do,]$min_ma/2),
     colSums(Obs1>0),type="o",xlim=rev(range(stages[do,]$max_ma)),xlab='Ma',ylab='# species',ylim=c(0,max(colSums(Obs1),colSums(Obs2))*1.02))
lines(rev(stages[do,]$max_ma/2+stages[do,]$min_ma/2),
      colSums(Obs2>0),type="o",xlim=rev(range(stages[do,]$max_ma)),col='red')

m1 <- make_BayesCMR(Obs1,dts,spec=~div,      ext=~d13C,   samp=~1,data=drivers)
m2 <- make_BayesCMR(Obs2,dts,spec=~d13C*div,   ext=~1,      samp=~SeaLev,data=drivers)
m12 <- make_BayesCMR_2clades(Obs1,Obs2,dts = dts,
                             spec1 = ~div1+div2+time,   ext1 = ~div1+div2+time,samp1 = ~time,
                             spec2 = ~div1+div2+time ,  ext2 = ~div1+div2+time,samp2 = ~time,data=drivers)

x = runif(m12$npar,min=-1.5,max=-1.4)
c(m1$probfun(x[seq(1,max(unlist(m12$inx$inx1)))]),
    m2$probfun(x[seq(1+max(unlist(m12$inx$inx1)),max(unlist(m12$inx$inx2)))]))
m12$probfun(x)

# Here is a bug in the Make_bayes for single clade; it makes the inx for samp drivers
# being wrongly indexes (and not counted). Think its fixed.

f1 <- MCMC_CMR(m1,niter=5e4)

f2 <- MCMC_CMR(m2,niter=1e5)
x0 = optim(runif(m12$npar,min=-0.01,max=0.1),m12$probfun,control=list(maxit=1e4,fnscale=-1))
f12 <- MCMC_CMR(m12,niter=1e5,draweps=5e2)

f12a <- contMCMC_CMR(f12)
# ,x0=c(runif(6,min=-1.8,max=-1.5),
#                          rep(0,m12$npar-6)),niter=2e4)


m12 <- make_BayesCMR_2clades(Obs1,Obs2,dts = dts,
                             spec1 = ~1,   ext1 = ~1,samp1 = ~1,
                             spec2 = ~1 ,  ext2 = ~1,samp2 = ~1,data=drivers)
f12 <- MCMC_CMR(m12,niter=1e4)


par(mfcol=c(2,2))
matplot(f1$Chain,type="l")
matplot(f2$Chain,type="l")
matplot(f12$Chain[,sort(unlist(m12$inx$inx1[1:3]))],type="l")
matplot(f12$Chain[,sort(unlist(m12$inx$inx2[1:3]))],type="l")



m1 <- make_BayesCMR(Obs1,dts,spec=~time+div,    ext=~time,   samp=~time  ,data=drivers)
m2 <- make_BayesCMR(Obs2,dts,spec=~time+div,    ext=~time,   samp=~time,data=drivers)
m12 <- make_BayesCMR_2clades(Obs1,Obs2,dts = dts,
                             spec1 = ~div12,   ext1 = ~1 ,  samp1 = ~time,
                             spec2 = ~div12,   ext2 = ~1,   samp2 = ~time,data=drivers)

x = runif(m12$npar,min=-1.5,max=-1.4)
sum(m1$probfun(x[seq(1,max(unlist(m12$inx$inx1)))]),
    m2$probfun(x[seq(1+max(unlist(m12$inx$inx1)),max(unlist(m12$inx$inx2)))]))
m12$probfun(x)


f12 <- MCMC_CMR(m12)
# Here is a bug in the Make_bayes for single clade; it makes the inx for samp drivers
# being wrongly indexes (and not counted). Fixed.
# It seems like things are working, when no ~time at least. need to test some more.

f1 <- MCMC_CMR(m1,niter=4e4)
f2 <- MCMC_CMR(m2,niter=4e4)
f12 <- MCMC_CMR(m12,x0=c(runif(6,min=-1.8,max=-1.5),
                         rep(0,m12$npar-6)),niter=2e5)
  par(mfcol=c(2,2))
matplot(f1$Chain[,1:3],type="l")
matplot(f2$Chain[,1:3],type="l")
matplot(f12$Chain[,sort(unlist(m12$inx$inx1[1:3]))],type="l")
matplot(f12$Chain[,sort(unlist(m12$inx$inx2[1:3]))],type="l")



# Seems legit

m12$likfun(x)$ll1$LogL
m1$likfun(x[1:3])$LogL
m12$likfun(x)$ll2$LogL
m2$likfun(x[4:6])$LogL

m12$lspecfun2(x)
m2$specfun(x[4:6])
m12$lextfun2(x)
m2$extfun(x[4:6])
m12$lsampfun2(x)
m2$sampfun(x[4:6])



# olderr <- options()$error
 options(error= browser())


