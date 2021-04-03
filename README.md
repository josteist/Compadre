
## Compadre - package to estimate drivers and diversification rates.

The Compadre package allows for the estimation of rates of speciation/origination, extinction and sampling using data from the fossil record. The approach is based on capture-mark-recapture techniques with additions to include potential drivers of rates of diversification, including diversity dependence. Definition of models occur by way of formulas and parameter estimation is Bayesian with internal functions to sample the model. Functions to plot rates and drivers are also included. The package also contains a selection of stage-level potential macroevolutionary drivers and a large dataset of Phanerozoic invertebrates downloaded from Paleobiology database.

## Installation
Download and install package from github, requires devtools
```
install.packages("devtools")
devtools::install_github("josteist/Compadre")
library(Compadre)
```

Compadre requires a couple of other packages (coda, mvtnorm etc), but you'll figure it out.

## Using Compadre

The basic input to a model is a matrix of 0/1 for observed taxa over intervals and a vector of temporal durations of these intervals. For each rate (speciation, extinction, and sampling) one can feed the model-generating function with basic formulas in R. Fossil matrices can easily be generated from [Paleobiology database](https://paleobiodb.org/), [NOW](http://www.helsinki.fi/science/now/ "NOW database") or other sources. 

First a simple example on diversification during the Jurassic & Cretaceous.

Extract the geological stages from the in-package timescale (`GSA_timescale`) and invertebrate fossil database (`Occ_species` and `Occ_genera`). Drivers are extracted from data.frame `Proxies`. `Occ_species` is the big observation matrix of species level data from Paleobiology database included in the package. A similar matrix (`Occ_genera`) is also attached for genera. See details in main reference.

```
# GSA_timescale is internal to the package. scale_level 5 is stages. Identical to timescale used by PBDB
stages = GSA_timescale[GSA_timescale$scale_level==5,] 

# Which intervals to use. 100 is Fortunian (first cambrian stage) 1 is Holocene
# Doing intervals indexed 45 through 23, Jurassic/Cretaceous
do = seq(45,23,by=-1) 

# Picking out the drivers for this example. All drivers are in object 'Proxies' internal to the package.
drivers <- Proxies[do,] 

# Getting the durations of the intervals
dts = stages[do,]$max_ma-stages[do,]$min_ma 

# Extracting the species level invert data to relevant times. 
Obs = (1*(Occ_species>0))[,do] 
# Removing the ones not observed in the Mesozoic
Obs <- Obs[which(rowSums(Obs)>0),] 
```

There are 15 021 species in this dataset
```
head(Obs) # Just to show the structure of the input.

                                      Hettangian Sinemurian Pliensbachian Toarcian Aalenian Bajocian Bathonian Callovian Oxfordian
Aviculopecten occidentalis                     0          0             0        0        0        0         0         0         0
Pleuromya uniformis                            0          1             1        1        1        1         1         1         1
Rouillieria michalkowii                        0          0             0        0        0        0         0         0         0
Pinna (Pinna) lanceolata                       0          0             0        0        0        0         0         0         1
Grammatodon (Cosmetodon) keyserlingii          0          0             0        0        0        0         0         0         1
Goniomya literata                              0          0             0        0        1        0         1         1         1


# Making a mid-point time reference for plotting
stages$midp = (stages$max_ma + stages$min_ma)/2
plot(stages[do,]$midp,colSums(Obs>0),type="o",ylab='Number of observed species',xlab='Million years ago',
     xlim=(rev(range(max(stages[do,]$max_ma)-cumsum(dts)))))
colbottom(stages[do,]) 
# This simple function adds colors at the bottom, from the stages data.frame

```

![raw_count](https://github.com/josteist/Compadre/blob/master/extra/fig1.png)

### Simple model
A simple model with no temporal variability and no drivers. Essentially: what is the mean speciation, extinction and sampling rates for Triassic/Jurassic? Any model is constructed by passing an obvervation matrix (taxa in rows, intervals in columns) with 1 indicating the taxa was observed within each interval. Additionally we pass a vector of durations for the intervals (here dts). After a model is generated the function MCMC_CMR can be used to sample the Bayesian posterior. 

```
# For reproducability.
set.seed('190480')
# Step 1: Make the model. Default is no drivers and no temporal variability in rates.
ma <- make_BayesCMR(Obs,dts)
# Step 2: Sample the model (It's bayesian, default priors)
fa <- MCMC_CMR(ma,niter=1e4)
# Written output:
summary(fa)
	 	 	=== Speciation rate parameters ===
              mean  2.5% median 97.5% p >/<0  Eff SS
Overall rate 0.107 0.105  0.107 0.109       0    225
	 	 	=== Extinction rate parameters ===
              mean  2.5% median 97.5% p >/<0  Eff SS
Overall rate 0.106 0.104  0.106 0.108       0    252
	 	 	=== Sampling rate parameters ===
              mean  2.5% median 97.5% p >/<0  Eff SS
Overall rate 0.151 0.146  0.151 0.155       0    316
```
Not surprisingly the rates of extinction and speciation are very similar. However, mean rates over many transitions are of course potentially very misleading measures of diversification. The last two columns show the probability of being 'significant', i.e. fraction of samples > or < 0, and effective sample sizes (says something about the convergence of the MCMC chain).

#### Plotting
Plot the rates. For models with no temporal variability and no drivers,
plots are just distributions of estimated mean rates.
```
plot(fa)
```

![plorate 1](https://github.com/josteist/Compadre/blob/master/extra/fig2.png)

### More complex model. 
Each rate (speciation, extinction and sampling) can also be modelled as a function of drivers. A driver can either be an external time-series of putative influences (temperature, sea level etc) or diversity dependence.

When generating the model using make_BayesCMR each rate can be defined by a formula `spec ~` ,`ext ~` and `samp ~`. 

` ~ time` (add a random deviation from the overall rate to each transition)
` ~ div`  (add a diversity-dependent term. 

Diversity is  estimated as # obs / sampling prob. Diversity can not impact sampling rates which are used for diversity estimation. An intercept for each rate is included by default. Also note that all drivers, including diversity, are normalized so parameters reflect relative importance as drivers. Lastly, all rates are implemented and estimated on the log-scale.

More complex models take longer to sample. Here speciation rates vary over time (they are not 'functions' of time as in + alpha*time, but this makes the model include a random effect per interval transition.)

```
mb <- make_BayesCMR(Obs,dts,
                    spec = ~time)
# Also you can call the model, and it will print out some basic info about it.
mb

  == Compadre model == 
Model includes  15021  taxa, spanning  23  intervals. 
The model was generated Sat Apr 03 08:09:43 2021 and has 26 parameters.
Model is of type  III . See ?make_BayesCMR for details.
 ==  Model terms  == 
Speciation ~ time 
Extinction ~ 1 
Sampling ~ 1 

# Fitting this model. 
fb <- MCMC_CMR(mb,niter=1e5) 
# niter gives the number of MCMC iterations. By defualt first half is burnin and it tunes a proposal distribution to generate good mixing. 
summary(fb)
	 	 	=== Speciation rate parameters ===
               mean   2.5% median 97.5% p >/<0  Eff SS
Overall rate 0.0993 0.0843 0.0985 0.118       0   62.7

	 	 	=== Extinction rate parameters ===
              mean 2.5% median 97.5% p >/<0  Eff SS
Overall rate 0.102  0.1  0.102 0.104       0    258

	 	 	=== Sampling rate parameters ===
              mean  2.5% median 97.5% p >/<0  Eff SS
Overall rate 0.144 0.139  0.144 0.149       0    231

plot(fb)
```

![fig fitfb](https://github.com/josteist/Compadre/blob/master/extra/fig3.png)

Effective sample size is a rough measure of whether or not the chains have converged,
i.e. if the parameters is estimated properly. A crude rule of thumb is
that the ESS should be > 200 for all parameters. ESS can be glanced from a summary for the main
parameters, and evaluated for all (including the random effects for each interval) by using the function ESS(fit).

Also note that the plot function outputs the sampled rates (and drawing of the plot is skipped by `drawplot=FALSE`)
```
tmp <- plot(fb,drawplot=FALSE)
names(tmp)
"SpecRates" "ExtRates"  "SampRates"

boxplot(t(tmp$SpecRates))
```
![figboxp](https://github.com/josteist/Compadre/blob/master/extra/fig4.png)

And you can also have normal y-axes for the internal rate-plots by setting `log=F` (`TRUE` by default). 
Additionally the function can be fed with a data.frame of the intervals used, from which it
extracts `$max_ma`, `$min_ma` (max and min age of interval) and `$color`. The `GSA_timescale` includes this
```
plot(fb,log=F,stages=stages[do,])
```
![figcolors](https://github.com/josteist/Compadre/blob/master/extra/fig5.png)

Here you see that the plots show one of the key assumptions of a CMR model; sampling rates are within intervals (plotted inside intervals), but speciation/extinction rates are for transitions (plotted ON the vertical lines separating intervals). Also, it is evident that doing analysis with temporally variable rate are of course to be preferred, given the data is sufficient. 

### A model including other drivers and diversity dependence
Each rate is defined as a separate function/formula and can also have interactions.

```
mc <- make_BayesCMR(Obs,dts,
                    spec = ~ div + d13C + time,
                    ext  = ~ d18O + time,
                    samp = ~ time,
                    data = drivers)
```
The drivers are extracted from the `data` inputted, but `div` and `time` are internal to the package and no entry in the data is needed (and also, probably best to avoid using variables with these names). All drivers a normalized to a mean of 0 and sd of 1 before estimation, so their relative importance is directly summarized by the estimated parameters. You could of course back-transform if you're interested in the actual effect sizes. Diversity is also normalized. All rates are estimated on the log-scale, so their impact is proportional, i.e. a parameter estimate of 1 leads to a doubling of speciation rate when the driver is 1 sd above mean, and a halfing of the rate when it's 1 sd below the mean.

Formulas are generic and also allows for interactions between the drivers. If you think that impact of diversity on speciation rates are more important with higher d13C values, then `spec = ~div*d13C` will estimate impacts of both drivers and their interaction. 

```
fc <- MCMC_CMR(mc,niter=1e6)
summary(fc)
	 	 	=== Speciation rate parameters ===
               mean    2.5% median  97.5% p >/<0  Eff SS
Overall rate  0.107  0.0841  0.106 0.1314  0.0000  142.6
div          -0.229 -0.4899 -0.237 0.0655  0.0587   60.6
d13C          0.194 -0.0590  0.198 0.4544  0.0713  118.7

	 	 	=== Extinction rate parameters ===
                mean    2.5%  median    97.5% p >/<0  Eff SS
Overall rate  0.0981  0.0709  0.0969  0.12755  0.0000    101
d18O         -0.2477 -0.5121 -0.2463 -0.00814  0.0233    128

	 	 	=== Sampling rate parameters ===
              mean  2.5% median 97.5% p >/<0  Eff SS
Overall rate 0.145 0.106  0.143 0.189       0    112

plot(fc,stages=stages[do,]) 

```

![figlastmod](https://github.com/josteist/Compadre/blob/master/extra/fig6.png)

The basic plot function outputs the estimated rates. plotDrivers outputs plots of the drivers, one figure for each rate.
```
plotDrivers(fc)
```

![figdriv1](https://github.com/josteist/Compadre/blob/master/extra/drivers1.png)
![figdriv2](https://github.com/josteist/Compadre/blob/master/extra/drivers2.png)

Here it seems like diversity has a negative impact on speciation rate and `d13C` a positive one, but the 95% credibility interval for both includes 0. There is a negative impact of `d18O` on extinctions rat, and in this case the credibility interval doees not include 0. Note that the effective sample size for some of the drivers are somewhat low, so here I would re-run with more iterations. 

#### Inspecting the MCMC-chain

You can also plot the chains to check for convergence. Large models might require relatively long runs, particularly with many drivers.
```
# First three are the 'mean' rates, on a log scale:
matplot(fc$Chain[,1:3],type="l",ylab='alpha - overall log rate',xlab='iteration')
```

![figchains](https://github.com/josteist/Compadre/blob/master/extra/chains.png)

If you want to know which parameter in the chain is which in the model then
```
mc$inx # prints their indexes into the chain.
$specInx
(Intercept)         div        d13C 
          1           4           5 

$extInx
(Intercept)        d18O 
          2           6 

$sampInx
(Intercept) 
          3 

$varInx
[1] 7 8 9

$specReInx
 [1] 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31

$extReInx
 [1] 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53

$sampReInx
 [1] 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74

```
So the samples for the `d18O` impact on extinction rate is in `fc$Chain[,6]`. `varInx` are indexes for the standard deviations of the random effects for speciation, extinction and sampling (here 7,8 and 9). `**ReInx` are the indexes for the random effects.


## Model selection
The approach also allows for model comparison by using Bayes factors. Using the same data we can compare a model with diversity-dependence (`md1`) vs one without (`md0`).

```
md0 <- make_BayesCMR(Obs,dts,
                    spec = ~ time,
                    ext  = ~ time,
                    samp = ~ time,
                    data = drivers)
md1 <- make_BayesCMR(Obs,dts,
                    spec = ~ div + time,
                    ext  = ~ time,
                    samp = ~ time,
                    data = drivers)
                    
fd0 <- MCMC_CMR(md0,niter=1e6)
fd1 <- MCMC_CMR(md1,niter=1e6)
```
To have sufficient evidence for the impact of diversity on speciation rates then i) the posterior distribution of the parameter should be different from 0 and ii) a model incorporating the driver should perform better in terms of model probability or likelihood. The package includes functions to estimate the Bayesian model likelihood by using importance sampling techniques. 

To get the log Bayesian Model Likelihood (BML) we use the function `doBML(CMR_fit)`. Since the BML is a sample from a distribution, a relatively large number of samples are needed. Here we draw 10 replicates for each model fit using default settings. The `replicate` below calls the function `doBML` 10 times and keeps the output `$logBML` for each replicate.

```
bml0 <- replicate(10,doBML(fd0)$logBML)
bml1 <- replicate(10,doBML(fd1)$logBML)

```
If we treat these two models as *a priori* equally likely, the Bayes factor when comparing them will be equal to the ratio of their Bayesian Model Likelihoods, and hence the different in their *log* BML. A different in logBML of 1-3 can be considered positive evidence, 3-5 strong evidence and a difference in logBML of > 5 can be seen as very strong (see Kass & Raftery 1995). A logBML difference of 3 corresponds to a Bayes Factor of approximately 20.

Plotting the draws of the logBML for both models above shows lack of evidence for the more complex model (including diversity-dependence). In this case there is insufficient evidence to conclude that diversity-dependence in rates of speciation is an important driver.

```
boxplot(data.frame(H0=bml0,DivDep=bml1),ylab='log Bayesian Model Likelihood')
grid()
```
![figbmls](https://github.com/josteist/Compadre/blob/master/extra/bmls.png)

## About the package
The general approach and methodological details are presented in **Compadre - estimating the drivers and dynamics of macroevolutionary change** submitted to Systematic Biology.

If you discover bugs or have wishes for improvements or changes, please do not hesitate to contact me. I'm also always in the mood for interesting collaborations.





