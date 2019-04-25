rm(list=ls())
library(RCurl)
library(Compadre)
stages <- GSA_timescale[GSA_timescale$scale_level==5,]
# THis URL call to Palebiology Databaser dataservice will download
# fossil observations for Cetacea - whales.
# Only finds resolved to species level and keep only accepted names.
js2 <- getURL('https://paleobiodb.org/data1.2/occs/list.csv?limit=all&base_name=Cetacea&idreso=species&pres=regular&show=acconly,phylo', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
Data_Iv <- read.csv(textConnection(js2),header=T) # Parsing the text
substr(js2,1,300)
# Download 080318: Eumetazoa^Vertebrata to species level
# 385095 obs.
# 250918 : 394903 obs
tmp = cbind(Data_Iv$max_ma,Data_Iv$min_ma) # span of observation
tmpst = c(stages$max_ma,1001) # boundaries
spans <- sapply(1:dim(Data_Iv)[1],function(i) min(which(tmp[i,1]<=tmpst))-min(which(tmp[i,2]<=tmpst)))
Data = Data_Iv[spans<3,];
Taxonomy<-(unique(Data[,c(14,15,16,17,18,8)]))
dim(Data) #  We have 1929 fossil observations
Taxonomy<-(unique(Data[,c(14,15,16,17,18,8)]))  # collating the Taxonomy
sapply(1:dim(Taxonomy)[2],function(ii){length(unique(Taxonomy[,ii]))})
# We have 1 phyla, 1 classes, 1 orders, 44 families,  351 genera and 590 species.

# All observations are put in a matrix with #species by #intervals.
#
unqspec <- unique(Data$accepted_no)
Oca <- t(sapply(1:length(unqspec),
              function(ii){hist(runif(sum(Data$accepted_no==unqspec[ii]),
                                      Data[Data$accepted_no==unqspec[ii],]$min_ma,
                                      Data[Data$accepted_no==unqspec[ii],]$max_ma),c(0,stages$max_ma,1000),
                                plot=F)$counts}))
# Keep only intervals with observations.
Oca = Oca[,1:100]
colnames(Oca) <- stages$interval_name
rownames(Oca) <- unique(Data[,c(6,8)])[,1]
Ocs = Oca[,rev(which(colSums(Oca)>0))]; # keep only intervals with obs
Ocs = Ocs[rowSums(Ocs)>0,]; #keep only species with obs
dts = rev((stages$max_ma-stages$min_ma)[which(colSums(Oca)>0)])

plot(sum(dts)-cumsum(dts),colSums(Ocs>0),type="o",xlim=c(50,0))


m1 <- make_BayesCMR(1*(Ocs>0),dts)
m2 <- make_BayesCMRv2(Ocs,dts)
f1 <- MCMC_CMR(m1,niter=2e4)
f2 <- MCMC_CMR(m2,niter=2e4)
plot(f1)
plot(f2)

m1a <- make_BayesCMR(Ocs>0,dts,RE=c(T,T,T))
m2a <- make_BayesCMRv2(Ocs>0,dts,RE=c(T,T,T))
f1a <- MCMC_CMR(m1a,niter=5e5,draweps=1e4)
f2a <- MCMC_CMR(m2a,niter=5e5,draweps=1e4)

plot(f1a)
plot(f2a)

boxplot(log10(t(sapply(seq(dim(f1a$Chain)[1]/2,dim(f1a$Chain)[1],by=10),
       function(ii){m1a$n_est(f1a$Chain[ii,])}))))
boxplot(log10(t(sapply(seq(dim(f1a$Chain)[1]/2,dim(f2a$Chain)[1],by=10),
                       function(ii){m2a$n_est(f2a$Chain[ii,])}))))
