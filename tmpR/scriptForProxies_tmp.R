#
# Bin averages.
rm(list=ls())

# Rationale is that d18O deviation from the equation proposed by veizer/Prokoph in 2015 form ain inverse proxy for temperature and 
# d13C isotopes signal productivity in that increased levels of d13C indicate higher productivity in the oceans (more sequestration of )

# mmc2 is the 'surface' data, which we will use.
# None of these are 'adjusted' in sheet 2.
vp2 <- as.data.frame(read_excel("C:/Users/josteist/Documents/Databases/VeizerProkoph2015/tmpmmc2.xlsx",sheet=2))
tmp2 <- data.frame(gts2004=as.numeric(vp2$gts2004),
                   gts2012=as.numeric(vp2$gts2012),
                   d13C=as.numeric(vp2$d13C),
                   d18O=as.numeric(vp2$d18O),
                   taxa=vp2$fossil,
                   clima = vp2$climate,
                   sht = rep(2,dim(vp2)[1]))
for (jj in 4:(length(excel_sheets("C:/Users/josteist/Documents/Databases/VeizerProkoph2015/tmpmmc2.xlsx")))){
  vp2 <- as.data.frame(read_excel("C:/Users/josteist/Documents/Databases/VeizerProkoph2015/tmpmmc2.xlsx",sheet=jj))
  tmp <- data.frame(gts2004=as.numeric(vp2$gts2004),
                     gts2012=as.numeric(vp2$gts2012),
                     d13C=as.numeric(vp2$d13C),
                     d18O=as.numeric(vp2$d18O),
                     taxa=vp2$fossil,
                     clima = vp2$climate,
                     sht = rep(jj,dim(vp2)[1]))
  tmp2 <- rbind(tmp2,tmp)
}
rm(tmp)
rm(vp2)

tmp2$taxa[which(tmp2$taxa %in% c('brachiopod','Brachiopod','Brachiopods'))] = 'Brachiopod'
# Remove all that are not Brachiopods of planktic Forams
tmp2<- tmp2[which(tmp2$taxa %in% c('Brachiopod','PlankticF')),]

# xyplot(as.numeric(d18)~gts2004 | taxa,data=tmp2)

# 1 - Locate all tropical and subtropical datapoints.
Alltropix <- unique(sort(c(grep('trop',unique(tmp2$clima),ignore.case=T),
                           grep('equ',unique(tmp2$clima),ignore.case=T),
                           grep('Deep',unique(tmp2$clima),ignore.case=T),
                           grep('Indik',unique(tmp2$clima),ignore.case=T))))
# 2. Which of these are subtropical and which are tropical.
subtropix <- intersect(unique(sort(c(grep('strop',unique(tmp2$clima),ignore.case=T),
                                     grep('subt',unique(tmp2$clima),ignore.case=T),
                                     grep('-s',unique(tmp2$clima),ignore.case=T),
                                     grep('trops',unique(tmp2$clima),ignore.case=T)))),
                       Alltropix)

tropix <- setdiff(Alltropix,subtropix)
# 3 - Find Artic and antarctic.
arctix <- unique(sort(c(grep('arct',unique(tmp2$clima),ignore.case=T),
                        grep('acti',unique(tmp2$clima),ignore.case=T))))
# 4 - Find temperate.
tempix <- unique(sort(c(grep('temp',unique(tmp2$clima),ignore.case=T),
                        grep('mid',unique(tmp2$clima),ignore.case=T),
                        grep('medit',unique(tmp2$clima),ignore.case=T),
                        grep('Atl-South',unique(tmp2$clima),ignore.case=T))))
# Are any sets found in more than 1 category?
diff(sort(c(subtropix,tropix,arctix,tempix)))
# Temp-Indik, temperate data from the Indian ocean I presume. Set to temperate

tropix = setdiff(tropix,9)
unique(tmp2$clima)[9]

tmp2$Climate = NA;
tmp2$Climate[which(tmp2$clima %in% unique(tmp2$clima)[tropix])] = 'Tropical'
tmp2$Climate[which(tmp2$clima %in% unique(tmp2$clima)[tempix])] = 'Temperate'
tmp2$Climate[which(tmp2$clima %in% unique(tmp2$clima)[arctix])] = 'Arctic/Antarctic'
tmp2$Climate[which(tmp2$clima %in% unique(tmp2$clima)[subtropix])] = 'Sub-tropical'

dataVP <- tmp2;
dataVP$taxa <- factor(dataVP$taxa)
dataVP$gts2012 <- as.numeric(dataVP$gts2012)
xyplot(d13C~ gts2012 | Climate,data=dataVP,group=taxa,xlim=c(540,-10),
       par.settings = list(superpose.symbol = list(pch = 19, cex = 0.8,
                                                   col = c(rgb(0.1,0.8,0,0.1),
                                                           rgb(0.8,0.2,0,0.1)))))
xyplot(d18O~ gts2012 | Climate,data=dataVP,group=taxa,xlim=c(540,-10),
       par.settings = list(superpose.symbol = list(pch = 19, cex = 0.8,
                                                   col = c(rgb(0.1,0.8,0,0.1),
                                                           rgb(0.8,0.2,0,0.1)))))

library(Compadre)
stages = GSA_timescale[GSA_timescale$scale_level==5,]; # Stages.

tmpdata <- subset(dataVP,Climate %in% c('Tropical','Sub-tropical'))
o_cor <- function(t){-0.00003*t^2 + 0.0046*t}
tmpdata$d18O = o_cor(tmpdata$gts2012)-tmpdata$d18O
bnavs <- sapply(1:dim(stages)[1],
                function(ii){
                  colMeans(tmpdata[which(tmpdata$gts2012<stages[ii,]$max_ma &
                                           tmpdata$gts2012>stages[ii,]$min_ma),3:4],na.rm=T)})

matplot(t(bnavs),type="l")
par(mfrow=c(1,2))
gm1 <- gam(d13C~s(gts2012),data=tmpdata)
gm2 <- gam(d13C~s(gts2012,bs="cr"),data=tmpdata)
gm3 <- gam(d13C~s(gts2012,k=70),data=tmpdata)
gm4 <- gam(d13C~s(gts2012,k=700),data=tmpdata)

gm1o <- gam(d18O~s(gts2012),data=tmpdata)
gm2o <- gam(d18O~s(gts2012,bs="cr"),data=tmpdata)
gm3o <- gam(d18O~s(gts2012,k=70),data=tmpdata)
gm4o <- gam(d18O~s(gts2012,k=700),data=tmpdata)

plot(tmpdata$gts2012,tmpdata$d13C,col=rgb(0,0,0.3,0.4),xlim=c(540,0),pch=20,xlab='Time',ylab=expression(paste(delta,"13C")))
lines(seq(0,540,by=1),predict(gm1,list(gts2012=seq(0,540,by=1))),type="l",col=rgb(1,0,0,0.6),lwd=3)
lines(seq(0,540,by=1),predict(gm2,list(gts2012=seq(0,540,by=1))),type="l",col=rgb(0,1,0,0.6),lwd=3)
lines(seq(0,540,by=1),predict(gm3,list(gts2012=seq(0,540,by=1))),type="l",col=rgb(0,0,1,0.6),lwd=3)
lines(stages$min_ma/2+stages$max_ma/2,bnavs[1,],type="l",lwd=3,col='grey')

legend('topright',c(expression(paste(delta,'18O')),'Stage means','k=9','k=70','k=700'),
       col=c(rgb(0,0,0.3,0.4),'grey',rgb(1,0,0,0.6),rgb(0,1,0,0.6),rgb(0,0,1,0.6)),
       lwd=4,pch=c(20,rep(NA,4)),lty=c(NA,1,1,1,1))


plot(tmpdata$gts2012,tmpdata$d18O,col=rgb(0,0,0.3,0.4),xlim=c(540,0),pch=20,xlab='Time',ylab=expression(paste(delta,"18O-deviations")))
lines(seq(0,540,by=1),predict(gm1o,list(gts2012=seq(0,540,by=1))),type="l",col=rgb(1,0,0,0.6),lwd=3)
lines(seq(0,540,by=1),predict(gm2o,list(gts2012=seq(0,540,by=1))),type="l",col=rgb(0,1,0,0.6),lwd=3)
lines(seq(0,540,by=1),predict(gm3o,list(gts2012=seq(0,540,by=1))),type="l",col=rgb(0,0,1,0.6),lwd=3)
lines(stages$min_ma/2+stages$max_ma/2,bnavs[2,],type="l",lwd=3,col='grey')


d13C <- rbind(bnavs[1,],sapply(1:dim(stages)[1],
                function(ii){c(mean(predict(gm1,list(gts2012=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,length.out=41)))),
                               mean(predict(gm2,list(gts2012=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,length.out=41)))),
                               mean(predict(gm3,list(gts2012=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,length.out=41)))),
                               mean(predict(gm4,list(gts2012=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,length.out=41)))))}))
rowSums((d13C - rep(d13C[1,],each=5))^2,na.rm=T)

  d18O <- rbind(bnavs[2,],sapply(1:dim(stages)[1],
               function(ii){c(mean(predict(gm1o,list(gts2012=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,length.out=41)))),
                              mean(predict(gm2o,list(gts2012=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,length.out=41)))),
                              mean(predict(gm3o,list(gts2012=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,length.out=41)))),
                              mean(predict(gm4o,list(gts2012=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,length.out=41)))))}))
rowSums((d18O[,1:90] - rep(d18O[1,1:90],each=5))^2,na.rm=T)
unlist(lapply(list(gm1o,gm2o,gm3o,gm4o),AIC))
# So what's the squared difference between stage-level means and the gams?


# So we will replace the missing bnavs with the ones from the GAM with k=70 for both isotope series

OI <- 



gm4o <- gam(d18O~s(gts2012,k=700,),data=tmpdata)


# so use BOTH "s" with min_ma and "S" with max_ma, to get lines
# at both ends.
# type="S" starts with two first ys at x[1], "s", draws y[1] at x[1]..x[2]
plot(tmpdata$gts2012,tmpdata$d18O,col=rgb(0,0,0.3,0.4),
     xlim=c(540,0),pch=20,xlab='Time',
     ylab=expression(paste(delta^"18","O-deviations")))
matplot(stages$min_ma,t(d18O),type="s",add=T,lwd=3,col=c('blue','red','green','purple'))
matplot(stages$max_ma,t(d18O),type="S",add=T)
abline(v=stages$min_ma,col=rgb(0.3,0.3,0.3,0.3))

plot(tmpdata$gts2012,tmpdata$d13C,col=rgb(0,0,0.3,0.4),xlim=c(480,0),pch=20,xlab='Time',ylab=expression(paste(delta,"13C")))
matplot(stages$min_ma,t(d13C),type="s",add=T,lwd=3,
        col=c('blue','red','green','purple'),lty=1)


matplot(t(d13C),type="l")
matplot(stages$min_ma,t(d13C),type="S")
matplot(stages$min_ma,t(d18O),type="S",xlim=c(540,0))

matplot(stages$min_ma/2+stages$max_ma/2,t(d13C),type="l",lwd=3,col=c('blue','red','green','purple'))

# 
# 
# 
# # CUTS BELOW, don't think we need GTS2004 scale, seems that all the data we want have been scaled to GTS2012
# table(tmp2$Climate)
# # We will only use tropical and sub-tropical (low latitude data)
# ix <- c( which(tmp2$Climate=='Tropical'),
#          which(tmp2$Climate=='Sub-tropical'),
#          which(tmp2$Climate=='Temperate'))
# DataIsotop <- tmp2[ix,];
# 
# # Gradstein timescale 2004
# GTS2004 <- c(0,0.0115,             0.126,             0.781,             1.806,             2.588,             3.6,             5.332,             7.246,             11.608,             13.65,             15.97,             20.43,             23.03,             28.4,             33.9,             37.2,             40.4,             48.6,             55.8,             58.7,             61.7,             65.5,             70.6,             83.5,             85.8,             89.3,             93.5,             99.6,             112,             125,             130,             136.4,             140.2,             145.5,             150.8,             155.7,             161.2,             164.7,             167.7,             171.6,             175.6,             183,             189.6,             196.5,             199.6,             203.6,             216.5,             228,             237,             245,             249.7,             251,             253.8,             260.4,             265.8,             268,             270.6,             275.6,             284.4,             294.6,             299,             303.9,             306.5,             311.7,             318.1,             326.4,             345.3,             359.2,             374.5,             385.3,             391.8,             397.5,             407,             411.2,             416,             418.7,             421.3,             422.9,             426.2,             428.2,             436,             439,             443.7,             445.6,             455.8,             460.9,             468.1,             471.8,             478.6,             488.3)
# # Also the current timescale (GTS 2012) from pbdb
# stages = GSA_timescale[GSA_timescale$scale_level==5,]; # Stages.
# GTS2012 <- unique(stages[1:91,]$min_ma);
# # ONly use the first
# which(is.na(DataIsotop$gts2012))
# tmp_gsa2012 <- sapply(which(is.na(DataIsotop$gts2012)),
#                       function(ii){rescale_time(GTS2004,
#                                                 GTS2012,DataIsotop$gts2004[ii])})
# 
# DataIsotop[which(is.na(DataIsotop$gts2012)),]$gts2012 = tmp_gsa2012
# # Now all with missing times on Gradstein 2012 is entered.
# 
# plot(DataIsotop$gts2012,DataIsotop$d13C,xlim=c(540,0))
# 
# tmpdata <- data.frame(time=as.numeric(DataIsotop$gts2012),
#                       d13C=as.numeric(DataIsotop$d13C),
#                       d18O=as.numeric(DataIsotop$d18O))
# tmpdata <- tmpdata[-which(is.na(tmpdata$d13C)),]
# 
# 
# 
# 
# 
# bnavs <- sapply(1:dim(stages)[1],
#        function(ii){
#          colMeans(tmpdata[which(tmpdata$time<stages[ii,]$max_ma &
#                             tmpdata$time>stages[ii,]$min_ma),],na.rm=T)})
# 
# # General additive models.
# gm1 <- gam(d13C~s(time),data=tmpdata)
# gm2 <- gam(d13C~s(time,bs="cr"),data=tmpdata)
# gm3 <- gam(d13C~s(time,k=99),data=tmpdata)
# plot(tmpdata$time,tmpdata$d13C,col=rgb(0,0,0.3,0.4),xlim=c(540,0),pch=20,xlab='Time',ylab=expression(paste(delta,"13C")))
# lines(seq(0,540,by=1),predict(gm1,list(time=seq(0,540,by=1))),type="l",col=rgb(1,0,0,0.6),lwd=3)
# lines(seq(0,540,by=1),predict(gm2,list(time=seq(0,540,by=1))),type="l",col=rgb(0,1,0,0.6),lwd=3)
# lines(seq(0,540,by=1),predict(gm3,list(time=seq(0,540,by=1))),type="l",col=rgb(0,0,1,0.6),lwd=3)
# 
# usefit <- gm3;
# gmavs <- sapply(1:dim(stages)[1],
#                 function(ii){
#                   mean(predict(usefit,list(time=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,by=0.1))))})
#                 })
# 
# plot(tmpdata$time,tmpdata$d13C,col=rgb(0,0,0.3,0.4),xlim=c(540,0),pch=20,xlab='Time',ylab=expression(paste(delta,"13C")))
# lines((stages$min_ma+stages$max_ma)/2,bnavs[2,],type="l",col=rgb(1,0,0,0.8),lwd=3)
# lines((stages$min_ma+stages$max_ma)/2,gmavs,type="l",col=rgb(0,1,0,0.8),lwd=3)
# 
# 
# # General additive models Oxygen
# gm1o <- gam(d18O~s(time),data=tmpdata)
# gm2o <- gam(d18O~s(time,bs="cr"),data=tmpdata)
# gm3o <- gam(d18O~s(time,k=99),data=tmpdata)
# usefit <- gm3o;
# gmavsO <- sapply(1:dim(stages)[1],
#                 function(ii){
#                   mean(predict(usefit,list(time=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,by=0.1))))})
# 
# 
# 
# plot(tmpdata$time,tmpdata$d18O,col=rgb(0,0,0.3,0.4),xlim=c(540,0),pch=20,xlab='Time',ylab=expression(paste(delta,"18O")))
# lines(seq(0,540,by=1),predict(gm1o,list(time=seq(0,540,by=1))),type="l",col=rgb(1,0,0,0.9),lwd=3)
# lines(seq(0,540,by=1),predict(gm2o,list(time=seq(0,540,by=1))),type="l",col=rgb(0,1,0,0.9),lwd=3)
# lines(seq(0,540,by=1),predict(gm3o,list(time=seq(0,540,by=1))),type="l",col=rgb(0,0,1,0.9),lwd=3)
# lines((stages$min_ma+stages$max_ma)/2,bnavs[3,],type="l",col=rgb(1,0,0,0.8),lwd=3,lty=3)
# lines((stages$min_ma+stages$max_ma)/2,gmavsO,type="l",col=rgb(0,1,0,0.8),lwd=3,lty=3)
# 
# par(mfrow=c(2,2))
# plot(tmpdata$time,tmpdata$d18O,col=rgb(0,0,0.3,0.4),xlim=c(540,0),pch=20,xlab='Time',ylab=expression(paste(delta,"18O")))
# lines(stages$min_ma/2+stages$max_ma/2,bnavs[3,],type="l",col=rgb(0.9,0.1,0.1,0.8),lwd=3,lty=3)
# lines(stages$min_ma/2+stages$max_ma/2,gmavsO,type="o",col=rgb(0.1,0.9,0.1,0.8),lwd=3,lty=3)
# plot(bnavs[3,],gmavsO,xlab='Bin means',ylab='GAM-means')
# plot(tmpdata$time,tmpdata$d13C,col=rgb(0,0.3,0,0.4),xlim=c(540,0),pch=20,xlab='Time',ylab=expression(paste(delta,"13C")))
# lines(stages$min_ma/2+stages$max_ma/2,bnavs[2,],type="l",col=rgb(0.9,0.1,0.1,0.8),lwd=3,lty=3)
# lines(stages$min_ma/2+stages$max_ma/2,gmavs,type="o",col=rgb(0.1,0.9,0.1,0.8),lwd=3,lty=3)
# plot(bnavs[2,],gmavs,xlab='Bin means',ylab='GAM-means')
# 
# 
# # WIth corectino from Veizer
# eq2 <- function(t){-0.00003*(t^2) + 0.0046*t}
# tmpdata$d18O_c <- eq2(tmpdata$time) - tmpdata$d18O
# gm1o <- gam(d18O_c~s(time,k=9),data=tmpdata)
# gm2o <- gam(d18O_c~s(time,k=99),data=tmpdata)
# gm3o <- gam(d18O_c~s(time,k=999),data=tmpdata)
# 
# usefit <- gm3o;
# gmavsO <- cbind(sapply(1:dim(stages)[1],function(ii){mean(predict(gm2o,list(time=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,by=0.1))))}),
#                 sapply(1:dim(stages)[1],function(ii){mean(predict(gm3o,list(time=seq(stages[ii,]$min_ma,stages[ii,]$max_ma,by=0.1))))}))
# 
# 
# for (ii in 1:dim(stages)[1]){
#                 function(ii){
#                   colMeans(tmpdata[which(tmpdata$time<stages[ii,]$max_ma &
#                                            tmpdata$time>stages[ii,]$min_ma),],na.rm=T)})
# }
# 
# plot(tmpdata$time,eq2(tmpdata$time)-tmpdata$d18O,col=rgb(0,0,0.3,0.4),xlim=c(540,0),ylim=c(-10,10),pch=20,xlab='Time',ylab=expression(paste(delta,"18O")))
# matplot((stages$min_ma+stages$max_ma)/2,gmavsO,type="l",col=rgb(0.5,.5,0,0.8),lwd=5,add=T)
# lines(seq(540,0,by=-1),predict(gm1o,list(time=seq(540,0,by=-1))),type="l",col=rgb(1,0,0,0.8),lwd=5,lty=3)
# 
# 
# # CUTS
# sapply(1:dim(stages)[1],
#        function(ii){
#          colMeans(as.numeric(DataIsotop[which(DataIsotop$gts2012<stages[ii,]$max_ma &
#                                      DataIsotop$gts2012>stages[ii,]$min_ma),c(3:5)])),na.rm=T)})
