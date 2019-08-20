# Downloading taxonomic opinions (current) from PBDB:
rm(list=ls())
library(RCurl)
library(mvtnorm)

## ===== GETTING SPATIAL AND TEMPORAL DIVISIONS ===  ====
# Prelims. Reading in interval information. Downloaded from PBDB earlier.
intervals <- read.csv('C:/Users/josteist/Documents/Databases/FromPaleoDB/intervals_PBDB.txt')
stages = intervals[intervals$scale_level==5,]
# Eumetazoa^Vertebrata
# Trying Eumetazoa - Vertebrata
# js2 <- getURL('https://paleobiodb.org/data1.2/occs/list.csv?limit=all&base_name=Trilobita&envtype=marine&idreso=genus&pres=regular&show=acconly,classext,subgenus', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
# So the subgenus and classext is useful. Trying with Eumetazoa^Vertebrata
# 30.07.2019
js <- getURL('https://paleobiodb.org/data1.2/occs/list.csv?limit=all&base_name=Eumetazoa^Vertebrata&envtype=marine&idreso=genus&pres=regular&show=acconly,classext,subgenus', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)

Data_raw <- read.csv(textConnection(js),header=T)
dim(Data_raw)
# There are 782410 occurrances. How many are species/genera level?
table(Data_raw$accepted_rank)

# Only keeping observations that span at most 3 geological stages.
tmp = cbind(Data_raw$max_ma,Data_raw$min_ma)
tmpst = c(stages$max_ma,1001) # stage boundaries.
spans <- sapply(1:dim(Data_raw)[1],
                function(i) min(which(tmp[i,1]<=tmpst))-min(which(tmp[i,2]<=tmpst)))
table(spans)
#
spans
#   0      1      2      3      4      5      6      7      8      9     10     11     12 
# 77493 523802 125354  24753  14634   5936   3899   1936   1932    803   1125    225    316 
# 
sum(spans<3)/length(spans)
# About 7 prct of obs cross >2 stage boundaries in temporal resolution. 

# Keep only occurrences that cross less than 3 boundaries
Data = Data_raw[spans<3,];
# Removed ~5 %,
# 371211 left
# Data = Data_Iv;
# Inspecting Taxonomy
Taxonomy<-(unique(Data[,c(14,15,16,17,18,8)]))
# We remove primarily terrestrial groups that are found
# in marine environments.
dim(Data)
Data = Data[-which(Data$class %in% c('Arachnida','Chilopoda','Cirripedia','Insecta','Myriapoda')),]
dim(Data)
# About 6600 removed. 

names(Data)[c(14,16,18,20,22,6)]
# 367449 obs
# 359822 observations. 368432 (OLD)
# redo Taxonomy
Taxonomy<-(unique(Data[,c(14,16,18,20,22,6)]))
Taxonomy_no <- unique(Data[,c(15,17,19,21,23,8)])
sapply(1:dim(Taxonomy)[2],    function(ii){length(unique(Taxonomy[,ii]))})
sapply(1:dim(Taxonomy_no)[2], function(ii){length(unique(Taxonomy_no[,ii]))})
# 300719:  27 / 88 / 460 / 3475 / 27340 / 87538
# dim(Taxonomy) yields 87560, i.e. some have same name, but not same classification?
# Using the _no as taxonomic information yields
#          27/ 88 / 461  / 3477 / 25094 / 87609
# We will use the numbers as unique identifiers of taxa.

# So we want two matrices; one for genera (subgenus elevated) and one for 
# species.

# THis is a much quicker OccArr generator
unqspec <- unique(Data[which(Data$accepted_rank=='species'),]$accepted_no)
Occ_species <- sapply(1:length(unqspec),
              function(ii){hist(runif(sum(Data$accepted_no==unqspec[ii]),
                                      Data[Data$accepted_no==unqspec[ii],]$min_ma,
                                      Data[Data$accepted_no==unqspec[ii],]$max_ma),c(0,stages$max_ma,1000),
                                plot=F)$counts})
# length(unqspec)
Taxonomy_species          <- data.frame(acc_no = unqspec)
tmp <- sapply(1:length(unqspec),function(ii){which(Data$accepted_no==unqspec[ii])[1]})
Taxonomy_species$phylum   <- Data[tmp,]$phylum
Taxonomy_species$class    <- Data[tmp,]$class 
Taxonomy_species$order    <- Data[tmp,]$order
Taxonomy_species$family   <- Data[tmp,]$family
Taxonomy_species$genus    <- Data[tmp,]$genus
Taxonomy_species$acc_name <- Data[tmp,]$accepted_name

# Removing all deemed to be observed before Cambrian, i.e. in row 101 are removed
rmv <- which(colSums(Occ_species[1:100,])==0)
Occ_species = Occ_species[1:100,-rmv];
# also removing from Taxonomy data.frame
Taxonomy_species = Taxonomy_species[-rmv,]

colnames(Occ_species) <- Taxonomy_species$acc_name
rownames(Occ_species) <- stages$interval_name
# save(Occ_species,Taxonomy_species,file='OccSpec_InvertPhanero_300719.RData')


# ====  Doing the genera. ====

# Making a genera array. Elevate all subgenera to genera.
# Keep genus_no if subgenus_no is NA, if not keep subgenus_no
tmp_array <- ifelse(!is.na(Data$subgenus_no),
                    Data$subgenus_no,Data$genus_no)
length(unique(tmp_array))
# 26776 unique genera with subgenus as genera.       
unqgen <- unique(tmp_array)
Occ_genera <- sapply(1:length(unqgen),
                      function(ii){hist(runif(sum(tmp_array==unqgen[ii]),
                                              Data[tmp_array==unqgen[ii],]$min_ma,
                                              Data[tmp_array==unqgen[ii],]$max_ma),c(0,stages$max_ma,1000),
                                        plot=F)$counts})
Taxonomy_genera          <- data.frame(acc_no = unqgen)
tmp <- sapply(1:length(unqgen),function(ii){which(tmp_array==unqgen[ii])[1]})
Taxonomy_genera$phylum   <- Data[tmp,]$phylum
Taxonomy_genera$class    <- Data[tmp,]$class 
Taxonomy_genera$order    <- Data[tmp,]$order
Taxonomy_genera$family   <- Data[tmp,]$family
Taxonomy_genera$genus    <- Data[tmp,]$genus
Taxonomy_genera$acc_name <- Data[tmp,]$accepted_name

rmv <- which(colSums(Occ_genera[1:100,])==0)
Occ_genera = Occ_genera[1:100,-rmv];
# also removing from Taxonomy data.frame
Taxonomy_genera = Taxonomy_genera[-rmv,]

colnames(Occ_genera) <- Taxonomy_genera$acc_name
rownames(Occ_genera) <- stages$interval_name



# save(Occ_genera,Taxonomy_genera,file='OccGenus_InvertPhanero_300719.RData')
# save.image('AllInverts_uncompiled_300719.RData')


