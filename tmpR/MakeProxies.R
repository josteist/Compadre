
# Making this into a Rmarkdown file,
# but update R


# This could perhaps be a markdown document as SI to show how we get all the proxies?

rm(list=ls())
# Making a script that captures ALL time series of environmental proxies::
library(RCurl)
library(Compadre)
stages = GSA_timescale[GSA_timescale$scale_level==5,]
# PNAS_201702297/FinalData/DiscreteTimeSeries.csv
tmptxt <- getURL('https://raw.githubusercontent.com/UW-Macrostrat/PNAS_201702297/master/FinalData/ContinuousTimeSeries.csv', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
Fragm <- read.csv(textConnection(tmptxt),header=T) # Parsing the text

# THis assumes stages is defined;
# stages <- GSA_timescale[GSA_timescale$scale_level==5,]
myfa <- splinefun(Fragm[,1],Fragm[,2],method='fmm')
FR_ZF <- sapply(1:dim(stages)[1],function(ii){
  mean(myfa(seq(stages[ii,]$min_ma,stages[ii,]$max_ma,by=0.01)))})
# SO this is the fragmentation index average inside stages
plot(c(-10,-10),xlim=c(500,0),ylim=c(0.3,.6),ylab='Fragmentation index')
for (ii in 1:100){
  # rect(stages[ii,]$max_ma,par("usr")[3],stages[ii,]$min_ma,par("usr")[3],col=toString(stages[ii,]$color))
  rect(stages[ii,]$max_ma,par("usr")[3],stages[ii,]$min_ma,400,col=paste0(toString(stages[ii,]$color),50),border=NA)

}
lines(Fragm[,1],Fragm[,2],type="l")
lines(stages$min_ma/2 + stages$max_ma/2, FR_ZF,col='red',type="o")

# Fragmentation done.


# Sea level curve, using composite from Cardenas & Harries 2010 Nature Geoscience.
# Supplementary is downloaded and pasted into a csv file
# Source "https://media.nature.com/original/nature-assets/ngeo/journal/v3/n6/extref/ngeo869-s1.pdf",

SL <- read.csv('C:/Users/josteist/Documents/Databases/AbioticProxyseries/ForRpack/Cardenas_Harries_TableS3_Sealevel_5myr.csv')

# Making a simple splineinterpolation.
# COuld try others...
myfa <- splinefun(SL[,1],SL[,2],method='fmm')
plot(c(-10,-10),xlim=c(500,0),ylim=c(0,100),ylab='Sea level curve')
for (ii in 1:100){
  rect(stages[ii,]$max_ma,par("usr")[3],stages[ii,]$min_ma,0,col=toString(stages[ii,]$color))
  rect(stages[ii,]$max_ma,par("usr")[3],stages[ii,]$min_ma,400,col=paste0(toString(stages[ii,]$color),50),border=NA)

  }

lines(seq(500,0,by=-0.1),myfa(seq(500,0,by=-0.1)),type="l")
points(SL[,1],SL[,2],col='red')
SL_CH <- sapply(1:dim(stages)[1],function(ii){
  mean(myfa(seq(stages[ii,]$min_ma,stages[ii,]$max_ma,by=0.01)))})



# Loading in d13C and d18O from equatorial and temperate areas from
# Veizer and Prokoph 2015.
setwd('C:/Users/josteist/Documents/Databases/VeizerProkoph2015')
library(readxl)
dir()



# mmc2 is the 'surface' data, which we will use.
# None of these are 'adjusted' in sheet 2.
vp2 <- as.data.frame(read_excel("C:/Users/josteist/Documents/Databases/VeizerProkoph2015/tmpmmc2.xlsx",sheet=2))
tmp2 <- data.frame(vp2$gts2004,vp2$gts2012,vp2$d13C,vp2$d18O,
                   vp2$`87sr/86Sr`,vp2$fossil,vp2$Location,vp2$climate)
for (jj in 4:(length(excel_sheets("C:/Users/josteist/Documents/Databases/VeizerProkoph2015/tmpmmc2.xlsx")))){
  vp2 <- as.data.frame(read_excel("C:/Users/josteist/Documents/Databases/VeizerProkoph2015/tmpmmc2.xlsx",sheet=jj))
  tmp <- data.frame(vp2$gts2004,vp2$gts2012,vp2$d13C,vp2$d18O,
                     vp2$`87sr/86Sr`,vp2$fossil,vp2$Location,vp2$climate)
  tmp2 <- rbind(tmp2,tmp)
}




# 1 - Locate all tropical and subtropical datapoints.
Alltropix <- unique(sort(c(grep('trop',unique(tmp2$vp2.climate),ignore.case=T),
                           grep('equ',unique(tmp2$vp2.climate),ignore.case=T),
                           grep('Deep',unique(tmp2$vp2.climate),ignore.case=T),
                           grep('Indik',unique(tmp2$vp2.climate),ignore.case=T))))
# 2. Which of these are subtropical and which are tropical.
subtropix <- intersect(unique(sort(c(grep('strop',unique(tmp2$vp2.climate),ignore.case=T),
                                     grep('subt',unique(tmp2$vp2.climate),ignore.case=T),
                                     grep('-s',unique(tmp2$vp2.climate),ignore.case=T),
                                     grep('trops',unique(tmp2$vp2.climate),ignore.case=T)))),
                       Alltropix)
unique(tmp2$vp2.climate)[subtropix]
unique(tmp2$vp2.climate)[setdiff(tropix,subtropix)]
tropix <- setdiff(Alltropix,subtropix)
# 3 - Find Artic and antarctic.
arctix <- unique(sort(c(grep('arct',unique(tmp2$vp2.climate),ignore.case=T),
                        grep('acti',unique(tmp2$vp2.climate),ignore.case=T))))
# 4 - Find temperate.
tempix <- unique(sort(c(grep('temp',unique(tmp2$vp2.climate),ignore.case=T),
                        grep('mid',unique(tmp2$vp2.climate),ignore.case=T),
                        grep('medit',unique(tmp2$vp2.climate),ignore.case=T),
                        grep('Atl-South',unique(tmp2$vp2.climate),ignore.case=T))))
# Are any sets found in more than 1 category?
diff(sort(c(subtropix,tropix,arctix,tempix)))
# Temp-Indik, temperate data from the Indian ocean I presume. Set to temperate

tropix = setdiff(tropix,9)
unique(tmp2$vp2.climate)[9]

tmp2$Climate = NA;
tmp2$Climate[which(tmp2$vp2.climate %in% unique(tmp2$vp2.climate)[tropix])] = 'Tropical'
tmp2$Climate[which(tmp2$vp2.climate %in% unique(tmp2$vp2.climate)[tempix])] = 'Temperate'
tmp2$Climate[which(tmp2$vp2.climate %in% unique(tmp2$vp2.climate)[arctix])] = 'Arctic/Antarctic'
tmp2$Climate[which(tmp2$vp2.climate %in% unique(tmp2$vp2.climate)[subtropix])] = 'Sub-tropical'


table(tmp2$Climate)
# We will only use tropical and sub-tropical (low latitude data)
ix <- c( which(tmp2$Climate=='Tropical'),
        which(tmp2$Climate=='Sub-tropical'),
        which(tmp2$Climate=='Temperate'))
DataIsotop <- tmp2[ix,];

# Gradstein timescale 2004
GTS2004 <- c(0,0.0115,             0.126,             0.781,             1.806,             2.588,             3.6,             5.332,             7.246,             11.608,             13.65,             15.97,             20.43,             23.03,             28.4,             33.9,             37.2,             40.4,             48.6,             55.8,             58.7,             61.7,             65.5,             70.6,             83.5,             85.8,             89.3,             93.5,             99.6,             112,             125,             130,             136.4,             140.2,             145.5,             150.8,             155.7,             161.2,             164.7,             167.7,             171.6,             175.6,             183,             189.6,             196.5,             199.6,             203.6,             216.5,             228,             237,             245,             249.7,             251,             253.8,             260.4,             265.8,             268,             270.6,             275.6,             284.4,             294.6,             299,             303.9,             306.5,             311.7,             318.1,             326.4,             345.3,             359.2,             374.5,             385.3,             391.8,             397.5,             407,             411.2,             416,             418.7,             421.3,             422.9,             426.2,             428.2,             436,             439,             443.7,             445.6,             455.8,             460.9,             468.1,             471.8,             478.6,             488.3)
# Also the current timescale (GTS 2012) from pbdb
stages = GSA_timescale[GSA_timescale$scale_level==5,]; # Stages.
GTS2012 <- unique(stages[1:91,]$min_ma);
# ONly use the first
which(is.na(DataIsotop$vp2.gts2012))
tmp_gsa2012 <- sapply(which(is.na(DataIsotop$vp2.gts2012)),
                      function(ii){rescale_time(GTS2004,
                                                GTS2012,DataIsotop$vp2.gts2004[ii])})

DataIsotop[which(is.na(DataIsotop$vp2.gts2012)),]$vp2.gts2012 = tmp_gsa2012
# Now all with missing times on Gradstein 2012 is entered.

plot(DataIsotop$vp2.gts2012,DataIsotop$vp2.d13C,xlim=c(540,0))

tmpdata <- data.frame(time=as.numeric(DataIsotop$vp2.gts2012),
                      d13C=as.numeric(DataIsotop$vp2.d13C),
                      d18O=as.numeric(DataIsotop$vp2.d18O))
tmpdata <- tmpdata[-which(is.na(tmpdata$d13C)),]

x = tmpdata$time;
y = tmpdata$d13C;

# THink we'll try stages-level means, instead of Gams.

plot(x,y)
# hwf to select these? the General cross validation value goes
# down, but how low can you go, lizm?
# Is setting k = -1 it optimizes it seems
gm1 <- gam(y~s(x,k=9))
gm2 <- gam(y~s(x,k=90))
gm3 <- gam(y~s(x,k=900))
gm4 <- gam(y~s(x,k=-1))
plot(gm1,ylim=c(-8,10))
points(x,y,col=rgb(0.1,0.1,0.1,0.1))

pre1 <- predict.gam(gm1,newdata=list(x=0:540))
pre2 <- predict.gam(gm2,newdata=list(x=0:540))
pre3 <- predict.gam(gm3,newdata=list(x=0:540))
pre4 <- predict.gam(gm4,newdata=list(x=0:540))
plot(x,y,col=rgb(0.1,0.1,0.1,0.5),pch=16,xlim=c(540,0))
lines(0:540,pre1,type="l",col='red',lwd=2,lty=1)
lines(0:540,pre2,type="l",col='green',lwd=2,lty=1)
lines(0:540,pre3,type="l",col='blue',lwd=2,lty=1)
lines(0:540,pre4,type="l",col='yellow',lwd=2)
# Well I wouldn't say that the optimal


IsotopStage <- findInterval(DataIsotop$vp2.gts2012,stages$min_ma)
# WHich interval do these belong to?
d13C_VP <- sapply(1:max(IsotopStage),
       function(jj){mean(as.numeric(DataIsotop[which(IsotopStage==jj),]$vp2.d13C),na.rm=T)})

mp = stages$min_ma/2+stages$max_ma/2
plot(mp[1:98],d13C_VP,type="o",xlim=c(520,0))
points(DataIsotop$vp2.gts2012,DataIsotop$vp2.d13C,pch='.',col=rgb(0.1,0.1,0.1,0.1))

# THere are some (at least 1 stage) missing here, and the
# values do not look like the ones in my Helsinki talk (som minor diffs)
jj=2
# WHat thef did I do here, I suspect I used some methodology that Tom Ezard
# used for one of his papers sometime...

# using mgcv package for GAM fitting.

?gam
tmpdat <- data.frame(
  time = as.numeric(DataIsotop$vp2.gts2012),
  d13C = as.numeric(DataIsotop$vp2.d13C))
tmpdat <- tmpdat[-which(is.na(tmpdat$d13C)),]
tmpg <- gam(d13C~s(time,4),data=tmpdat,knots=10)




for (jj in 2:17){
  vp1 <- read_excel("1-s2.0-S0012825215000604-mmc2.xls",sheet=jj)
  print(names(vp1))
}
# So this sheet also has strontium, both also have d13C
# use 'adjusted?'


# To allow for import two changes to the raw files from Veizer& Prokoph 2015 was made:
# - both appendix 1 & 2 were renamed and changed into xlsx format, because it could not be importet for some odd reason.
# - appendix 2 lacked one some uniformity in names and columns, and only select columns are stored.



vp1 <- as.data.frame(read_excel("tmpmmc1.xlsx",sheet=2))
for (jj in 3:6){
  tmp <- as.data.frame(read_excel("tmpmmc1.xlsx",sheet=jj))
  tmp1 <- vp1
  vp1 <- rbind(tmp1,tmp)
}
# When adjusted values are not given, then the reported are used.
vp1$d18O_all = vp1$d18O_adj
vp1$d18O_all[is.na(vp1$d18O_all)] = vp1$d18O[is.na(vp1$d18O_all)]

vp1$d13C_all = vp1$d13C_adj
vp1$d13C_all[is.na(vp1$d13C_all)] = vp1$d13C[is.na(vp1$d13C_all)]


# None of these are 'adjusted' in sheet 2.
tmp2 <- as.data.frame(read_excel("tmpmmc2.xlsx",sheet=2))
vp2 <- data.frame(tmp2$gts2004,tmp2$gts2012,tmp2$d13C,tmp2$d18O,tmp2$`87sr/86Sr`,tmp2$fossil,tmp2$Location,tmp2$climate)
for (jj in 4:(length(excel_sheets("tmpmmc2.xlsx")))){
  tmp2 <- as.data.frame(read_excel("tmpmmc2.xlsx",sheet=jj))
  tmp <- data.frame(tmp2$gts2004,tmp2$gts2012,tmp2$d13C,tmp2$d18O,tmp2$`87sr/86Sr`,tmp2$fossil,tmp2$Location,tmp2$climate)
  vp2 <- rbind(vp2,tmp)
}

# THe climate zones are highly inconsistent, but there is a pattern.
# Collecting all that are "tropical" in some sense (some denoted equatorial in appendix from Verizer & Prokoph)
# Some of these are prob sub-tropical.


# 1 - Locate all tropical and subtropical datapoints.
Alltropix <- unique(sort(c(grep('trop',unique(tmp2$vp2.climate),ignore.case=T),
  grep('equ',unique(tmp2$vp2.climate),ignore.case=T),
  grep('Deep',unique(tmp2$vp2.climate),ignore.case=T),
  grep('Indik',unique(tmp2$vp2.climate),ignore.case=T))))
# 2. Which of these are subtropical and which are tropical.
subtropix <- intersect(unique(sort(c(grep('strop',unique(tmp2$vp2.climate),ignore.case=T),
                        grep('subt',unique(tmp2$vp2.climate),ignore.case=T),
                        grep('-s',unique(tmp2$vp2.climate),ignore.case=T),
                        grep('trops',unique(tmp2$vp2.climate),ignore.case=T)))),
                       Alltropix)
unique(tmp2$vp2.climate)[subtropix]
unique(tmp2$vp2.climate)[setdiff(tropix,subtropix)]
tropix <- setdiff(Alltropix,subtropix)
# 3 - Find Artic and antarctic.
arctix <- unique(sort(c(grep('arct',unique(tmp2$vp2.climate),ignore.case=T),
               grep('acti',unique(tmp2$vp2.climate),ignore.case=T))))
# 4 - Find temperate.
tempix <- unique(sort(c(grep('temp',unique(tmp2$vp2.climate),ignore.case=T),
                 grep('mid',unique(tmp2$vp2.climate),ignore.case=T),
                 grep('medit',unique(tmp2$vp2.climate),ignore.case=T),
                 grep('Atl-South',unique(tmp2$vp2.climate),ignore.case=T))))
# Are any sets found in more than 1 category?
diff(sort(c(subtropix,tropix,arctix,tempix)))
# Temp-Indik, temperate data from the Indian ocean I presume. Set to temperate

tropix = setdiff(tropix,9)
unique(tmp2$vp2.climate)[9]

tmp2$Climate = NA;
tmp2$Climate[which(tmp2$vp2.climate %in% unique(tmp2$vp2.climate)[tropix])] = 'Tropical'
tmp2$Climate[which(tmp2$vp2.climate %in% unique(tmp2$vp2.climate)[tempix])] = 'Temperate'
tmp2$Climate[which(tmp2$vp2.climate %in% unique(tmp2$vp2.climate)[arctix])] = 'Arctic/Antarctic'
tmp2$Climate[which(tmp2$vp2.climate %in% unique(tmp2$vp2.climate)[subtropix])] = 'Sub-tropical'

# All datapoints in tmp2  have been classified.

# vp2 = tmp2;


# Some points have not been rescaled to 2012 timescale.
# REscaling to updated Gradstein timescale.
# Importing old timescale for isotope data: GTS 2004. This is undefined some some of the cambrian boundaries, but for the data here that is missing GTS2012 this is not relevant
GTS2004 <- c(0,0.0115,             0.126,             0.781,             1.806,             2.588,             3.6,             5.332,             7.246,             11.608,             13.65,             15.97,             20.43,             23.03,             28.4,             33.9,             37.2,             40.4,             48.6,             55.8,             58.7,             61.7,             65.5,             70.6,             83.5,             85.8,             89.3,             93.5,             99.6,             112,             125,             130,             136.4,             140.2,             145.5,             150.8,             155.7,             161.2,             164.7,             167.7,             171.6,             175.6,             183,             189.6,             196.5,             199.6,             203.6,             216.5,             228,             237,             245,             249.7,             251,             253.8,             260.4,             265.8,             268,             270.6,             275.6,             284.4,             294.6,             299,             303.9,             306.5,             311.7,             318.1,             326.4,             345.3,             359.2,             374.5,             385.3,             391.8,             397.5,             407,             411.2,             416,             418.7,             421.3,             422.9,             426.2,             428.2,             436,             439,             443.7,             445.6,             455.8,             460.9,             468.1,             471.8,             478.6,             488.3)
# Also the current timescale (GTS 2012) from pbdb
stages = GSA_timescale[GSA_timescale$scale_level==5,]; # Stages.
GTS2012 <- unique(stages[1:91,]$min_ma);
# ONly use the first
which(is.na(vp1$`Age (GTS2012)`))

tmp_gsa2014 <- sapply(which(is.na(tmp2$vp2.gts2012)),
                      function(ii){rescale.time(GTS2004,
                                                GTS2012,tmp2$vp2.gts2004[ii])})
tmp2[which(is.na(tmp2$vp2.gts2012)),]$vp2.gts2012 <- tmp_gsa2014






plot(tmp2$vp2.gts2012,tmp2$vp2.d18O,xlim=c(540,0),
     col=rgb(0.2,0.2,0.2,0.2))
points(vp1$`Age (GTS2012)`,vp1$d18O_all,col=rgb(0.2,0.2,0.4,0.2))

plt <- rainbow(12)
for (jj in 1:length(unique(vp1$`Ocean:Region`))){
  ixtmp <- which(vp1$`Ocean:Region`==unique(vp1$`Ocean:Region`)[jj])
  points(vp1$`Age (GTS2012)`[ixtmp],vp1$d18O_all[ixtmp],
         col=plt[jj])

}

# par(mfrow=c(2,2))
for (jj in 1:length(unique(tmp2$Climate))){
  ixtmp <- which(tmp2$Climate == unique(tmp2$Climate)[jj]);
  points(tmp2$vp2.gts2012[ixtmp],tmp2$vp2.d18O[ixtmp],
         col=plt[jj+5],xlim=c(540,0))
}


# par(mfrow=c(2,2))
plot(tmp2$vp2.gts2012,tmp2$vp2.d13C,
     col=rgb(0.3,0.3,0.3,0.2),xlim=c(540,0))

for (jj in 1:length(unique(tmp2$Climate))){
  ixtmp <- which(tmp2$Climate == unique(tmp2$Climate)[jj]);
  points(tmp2$vp2.gts2012[ixtmp],tmp2$vp2.d13C[ixtmp],
       col=plt[jj],xlim=c(540,0))
}


Proxies <- array(NA,c(dim(stages)[1],4))

## =?===== TO  DO
# Sort all the tmp2 data into tropical and subtropical + tropical
# and do interval means. THese will be our 'proxies'


# And, everything in vp1 is 'Deep sea' data, and we
# are most likely interested in 'Surface' data. Actually, i think
# we should mostly use the tropical (perhaps also subtrop) surface data.
# PLotting using tropical and tropical + subtropitcal
ixtmp2 <- findInterval(tmp2$vp2.gts2012,stages$min_ma)
ixtmp1 <- findInterval(vp1$`Age (GTS2012)`,stages$min_ma)
tmpd13C_T <- sapply(1:100,function(ii){
  mean(c(as.numeric(tmp2$vp2.d13C[which(ixtmp2==ii)]),
         as.numeric(vp1$d13C_all[which(ixtmp1==ii)])),na.rm=T)})

tmpd18O_T <- sapply(1:100,function(ii){
  mean(c(as.numeric(tmp2$vp2.d18O[which(ixtmp2==ii)]),
         as.numeric(vp1$d18O_all[which(ixtmp1==ii)])),na.rm=T)})

#pa
par(mfrow=c(2,2))
plot(vp1$`Age (GTS2004)`,vp1$d13C_all,xlim=c(500,0),col=rgb(0.2,0.5,0.5),ylim=c(-5,10))
points(tmp2$vp2.gts2004,tmp2$vp2.d13C)
lines(stages$max_ma/2+stages$min_ma/2,
      tmpd13C_T,type="l",col='grey',lty=1,lwd=3)


plot(vp1$`Age (GTS2004)`,vp1$d18O_all,xlim=c(540,0),col=rgb(0.2,0.2,0.2,0.2),ylim=c(-14,4))
points(tmp2$vp2.gts2004,tmp2$vp2.d18O)
lines(stages$max_ma/2+stages$min_ma/2,
      tmpd18O_T,type="l",col='grey',lty=1,lwd=3)


# TROPICAL + SUB-TROPICAL
ixtmp2 <- findInterval(tmp2$vp2.gts2012,stages$min_ma)
ixtmp1 <- findInterval(vp1$`Age (GTS2012)`,stages$min_ma)
tmpd13C_TST <- sapply(1:100,function(ii){
  mean(as.numeric(tmp2[which(ixtmp2==ii
                                        & c(tmp2$Climate=='Tropical' ||
                                              tmp2$Climate=='Sub-tropical')),]$vp2.d13C),na.rm=T)})
tmpd18O_TST <- sapply(1:100,function(ii){
  mean(as.numeric(tmp2[which(ixtmp2==ii & c(tmp2$Climate=='Tropical' ||
                                              tmp2$Climate=='Sub-tropical')),]$vp2.d18O),na.rm=T)})

# PLotting tropical and tropical + sub-tropical
par(mfrow=c(2,2))
plot(tmp2$vp2.gts2012[c(which(c(tmp2$Climate=='Tropical' ||
                                  tmp2$Climate=='Sub-tropical'))],
     tmp2$vp2.d18O[c(which(tmp2$Climate=='Tropical'),
                     which(tmp2$Climate=='Sub-tropical'))],xlim=c(540,0),col=rgb(0.7,0,0,0.1),
     xlab='Ma',ylab='d18O')
lines(stages$max_ma/2+stages$min_ma/2,tmpd18O_TST,type='o',lwd=2)

plot(tmp2$vp2.gts2012[which(tmp2$Climate=='Tropical')],
     tmp2$vp2.d13C[which(tmp2$Climate=='Tropical')],col=rgb(0,0.5,0,0.1),xlim=c(540,0),
     xlab='Ma',ylab='d13C')
lines(stages$max_ma/2+stages$min_ma/2,tmpd13C_TST,type='o',lwd=2)






tmpd13C <- sapply(1:100,function(ii){
  mean(as.numeric(tmp2[which(ixtmp2==ii
                             & tmp2$Climate=='Tropical'),]$vp2.d13C),na.rm=T)})
tmpd18O <- sapply(1:100,function(ii){
  mean(as.numeric(tmp2[which(ixtmp2==ii & tmp2$Climate=='Tropical'),]$vp2.d18O),na.rm=T)})


par(mfrow=c(2,1))
plot(tmp2$vp2.gts2012[which(tmp2$Climate=='Tropical')],
     tmp2$vp2.d18O[which(tmp2$Climate=='Tropical')],xlim=c(540,0),col=rgb(0.7,0,0,0.1),
     xlab='Ma',ylab='d18O')
lines(stages$max_ma/2+stages$min_ma/2,tmpd18O,type='o',lwd=2)

plot(tmp2$vp2.gts2012[which(tmp2$Climate=='Tropical')],
     tmp2$vp2.d13C[which(tmp2$Climate=='Tropical')],col=rgb(0,0.5,0,0.1),xlim=c(540,0),
     xlab='Ma',ylab='d13C')
lines(stages$max_ma/2+stages$min_ma/2,tmpd13C,type='o',lwd=2)


# So it seems like the vp1 data is only recent, right? Yes, from 110 mya
# or something.


writeClipboard(as.character(unique(tmp2$vp2.climate)))

xyplot(vp2.d13C~vp2.gts2012 | Climate,data=tmp2)

# vp1 has some adjusted and some non-adjusted values. If adjusted exists,
# we use these, if not we use the ones reported (older data is adjusted).
d13C <- data.frame(ma = c(vp1$`Age (GTS2004)`[is.na(vp1$d13C_adj)],
                          vp1$`Age (GTS2004)`[!is.na(vp1$d13C_adj)],
                          tmp2$vp2.gts2004),
                   d13C = c(vp1$d13C[is.na(vp1$d13C_adj)],
                            vp1$d13C_adj[!is.na(vp1$d13C_adj)],
                            tmp2$vp2.d13C))




plot(vp1$`Age (GTS2004)`,vp1$d18O,pch='.',col=rgb(0.9,.8,0.1),xlim=c(500,0),ylim=c(-15,5))
points(tmp2$vp2.gts2004,tmp2$vp2.d18O,,pch='.',col=rgb(0.9,.7,0.1))

plot(vp1$`Age (GTS2004)`,vp1$d13C,pch='.',col=rgb(0.1,.5,.5),xlim=c(500,0),ylim=c(-4,8),ylab='{delta}_1_3C')
points(tmp2$vp2.gts2004,tmp2$vp2.d13C,pch='.',col=rgb(0.1,0.8,0.4))

plot(tmp2$vp2.gts2004,tmp2$vp2..87sr.86Sr.,col=rgb(0.3,0.3,0),xlim=c(500,0))


unique(vp1$climate)





