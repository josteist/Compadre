
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

The basic input to a model is a matrix of 0/1 for observed taxa over intervals and a vector of temporal durations of these intervals. For each rate (speciation, extinction, and sampling) one can feed the model-generating function with basic formulas in R. Fossil matrices can easily be generates from [Paleobiology database](https://paleobiodb.org/), [NOW](http://www.helsinki.fi/science/now/ "NOW database") or other sources. 

First a simple example on diversification during the Jurassic & Cretaceous.

Extract the geological stages from the in-package timescale (`GSA_timescale`) and invertebrate fossil database (`InvertPBDB`). Drivers are extracted from data.frame `Proxies`.
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

# Extracting the species level invert data to relevant times. InvertPBDB is the big matrix included in the package.
Obs = (1*(InvertPBDB>0))[,do] 
# Removing the ones not observed in the Mesozoic
Obs <- Obs[which(rowSums(Obs)>0),] 
```

There are 14912 species in this dataset.
```
head(Obs) # Just to show the structure of the input.

                           Hettangian Sinemurian Pliensbachian Toarcian Aalenian Bajocian
Aviculopecten occidentalis          0          0             0        0        0        0
Pleuromya uniformis                 0          1             1        1        1        1
Rouillieria michalkowii             0          0             0        0        0        0
Pinna (Pinna) lanceolata            0          0             0        0        0        0
Goniomya literata                   0          0             0        1        0        0
Idonearca capax                     0          0             0        0        0        0


plot(max(stages[do,]$max_ma)-cumsum(dts),colSums(Obs>0),type="o",ylab='# species',xlab='Ma',
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
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				  0.18	  0.17	  0.18	  0.18	  0.00	407.97	

	 	 	=== Extinction rate parameters ===
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				  0.18	  0.17	  0.18	  0.18	  0.00	407.44	

	 	 	=== Sampling rate parameters ===
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				  0.15	  0.14	  0.15	  0.15	  0.00	414.58	
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

Diversity is roughly estimated as # obs / sampling prob). Doesn't work for sampling rates, which are used for diversity estimation. An intercept is included by default. Also note that all drivers, including diversity, are normalized so parameters reflect relative importance as drivers. Lastly, all rates are implemented and estimated on the log-scale

More complex models take longer to estimate. Here speciation rates vary over time (they are not 'functions' of time as in + alpha*time, but this makes
the model include a random effect per interval transition.)

```
mb <- make_BayesCMR(Obs,dts,
                    spec = ~time)
# Also you can call the model, and it will print out some basic info about it.
mb

 == Compadre model == 
Model includes  14912  taxa, spanning  23  intervals. 
The model was generated Tue Jul 09 15:17:10 2019 and has 26 parameters.

 ==  Model terms  == 
Speciation ~ time 
Extinction ~ 1 
Sampling ~ 1 

# Fitting this model. 
fb <- MCMC_CMR(mb,niter=1e5) 
# niter gives the number of MCMC iterations. By defualt first half is burnin and it tunes a proposal distribution to generate good mixing. 
summary(fb)
	 	 	=== Speciation rate parameters ===
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				  0.17	  0.14	  0.17	  0.19	  0.00	119.45	

	 	 	=== Extinction rate parameters ===
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				  0.17	  0.17	  0.17	  0.18	  0.00	182.02	

	 	 	=== Sampling rate parameters ===
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				  0.14	  0.14	  0.14	  0.15	  0.00	316.56	

plot(fb)
```

![fig fitfb](https://github.com/josteist/Compadre/blob/master/extra/fig3.png)

Effective sample size is a rough measure of whether or not the chains have converged,
i.e. if the parameters is estimated properly. A crude rule of thumb is
that the ESS should be > 200 for all parameters. ESS can be glanced from a summary for the main
parameters, and evaluated for all (including the random effects for each interval) by using the function ESS(fit).

Also note that the plot function outputs the sampled rates;
```
tmp <- plot(fb)
names(tmp)
"SpecRates" "ExtRates"  "SampRates"

boxplot(t(tmp$SpecRates))
```
![figboxp](https://github.com/josteist/Compadre/blob/master/extra/fig3.png)

And you can also have normal y-axes for the internal rate-plots by setting `log=F` (`TRUE` by default). 
Additionally the function can be fed with a data.frame of the intervals used, from which it
extracts `$max_ma`, `$min_ma` (max and min age of interval) and `$color`. The `GSA_timescale` includes this
```
plot(fb,log=F,stages=stages[do,])
```
![figcolors](https://github.com/josteist/Compadre/blob/master/extra/fig5.png)

Here you see that the plots show one of the key assumptions of a CMR model; sampling rates are within intervals, but speciation/extinction rates are for transitions, i.e. going from one interval to next (plotted ON the stage limits in grey). Also, it is evident that doing analysis with temporally variable rate are of course to be preferred, given the data is sufficient. 

### A model including other drivers and diversity dependence
Each rate is defined as a separate function/formula and can also have interactions.

```
mc <- make_BayesCMR(Obs,dts,
                    spec = ~ div + d13C + time,
                    ext  = ~ d180_cor + time,
                    samp = ~ time,
                    data = drivers)
```
The drivers are extracted from the data inputted, but `div` and `time` are internal to the package and no entry in the data is needed (and also, please
avoid using variables with these names). All drivers a normalized to a mean of 0 and sd of 1 before estimation, so their relative importance is directly summarized by the estimated parameters. You could of course pack-transform if you're interested in the actual effect sizes. Diversity is also normalized. All rates are estimated on the log-scale, so their impact is proportional, i.e. a parameter estimate of 1 leads to a doubling of speciation rate when the driver is 1 sd above mean, and a halfing of the rate when it's 1 sd below the mean.

Formulas are generic and also allows for interactions between the drivers. If you think that impact of diversity on speciation rates are more important 
with higher d13C values, then spec = ~div*d13C will estimate impacts of both drivers and their interaction. 

```
fc <- MCMC_CMR(mc,niter=1e6)
summary(fc)
	 	 	=== Speciation rate parameters ===
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				 0.14	 0.11	 0.14	 0.18	 0.00	103.17	
div					-0.45	-0.68	-0.45	-0.21	 0.00	38.37	
d13C					3.8e-01	1.1e-01	3.7e-01	7.3e-01	2.1e-03	1.3e+02	

	 	 	=== Extinction rate parameters ===
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				  0.17	  0.12	  0.17	  0.25	  0.00	119.81	
d180_cor				 -0.212	 -0.553	 -0.208	  0.108	  0.096	264.343	

	 	 	=== Sampling rate parameters ===
       					mean	2.5%	median	97.5%	p >/<0 	Eff SS	
Overall rate				  0.131	  0.097	  0.129	  0.174	  0.000	148.828	

plot(fc) 

```

![figlastmod](https://github.com/josteist/Compadre/blob/master/extra/fig6.png)

The basic plot function outputs the estimated rates. plotDrivers outputs plots of the drivers, one figure for each rate.
```
plotDrivers(fc)
```

![figdriv1](https://github.com/josteist/Compadre/blob/master/extra/drivers1.png)
![figdriv2](https://github.com/josteist/Compadre/blob/master/extra/drivers2.png)

Here it seems like diversity has a negative impact on speciation rate, `d13C` a positive one. And perhaps a positive impact of `d180_cor` on extinctions rate, but most likely not significant (or rather the credibility interval includes 0). Also note that the effective sample size for `d13C` is rather low, so here I would re-run with more iterations. 

#### Inspecting the MCMC-chain

You can also plot the chains to check for convergence. Large models might require relatively long runs, particularly with many drivers.
```
dim(fc$Chain)
# 100000 samples of 95 parameters
# First three are the 'mean' rates:
matplot(fc$Chain[,1:3],type="l")
```

![figchains](https://github.com/josteist/Compadre/blob/master/extra/chains.png)

If you want to know which parameter in the chain is which in the model then
```
mc$inx # prints their indexes into the chain.
$specInx
(Intercept)         div        d13C 
          1           4           5 

$extInx
(Intercept)    d180_cor 
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


```


## In development.
This package is not entirely finished yet, but most functionality should be there. If you experience troubles, have tips or comments, feel free to contact me. 

