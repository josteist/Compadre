#' Function to simulate a birth-death-fossilization process.
#'
#' This function simulates a birth-death process with a fossilization/sampling scheme. Rates of speciation and extinction
#' can be dependent on time \emph{t} and standing diversity/richness, \emph{n}. Sampling rates can be dependent om time.
#' @param spec Speciation rate. Either given as a fixed rate or as a function of time (\emph{t}) \emph{and} richness (\emph{n}). Defaults to 0.1.
#' @param ext Extinction rate. Given as fixed rate or function of time and richness. Defaults to 0.001.
#' @param samp Sampling rate, given as fixed rate or function of time (\emph{t}) only. Defaults to 0.3
#' @param n_init Initial number of lineages. Defaults to 100.
#' @param dt_ints Array of 'intervals' in which to generate the fossil record, as an array of interval durations. All lineages sampled within these intervals are placed in the interval, regardless of precise time it was sampled for the output \emph{FosRec}. Defaults to rep(1,10)
#' @return A named list with \emph{Taxa} (time of origin and extinction of all taxa), \emph{Foss} (list of all fossilization/sampling event for all taxa), \emph{FosRec} (occurrence matrix with dimensions sampled species by intervals). The remainding entries in the list are the inputs given to the function, \emph{dts, Spec, Ext, Samp, n_init}
#'
#' @export

sim_BD <- function(spec=0.1,
                   ext=0.02,
                   samp = 0.3,
                   n_init=100,dt_ints=rep(1,10)){
  # Simulating a birth death process with spec and ext
  # as functions of (time, n_species now). Output
  # is matrix of origination and extinction times. Speciation is budding.
  #
  if (any(sapply(list(spec,ext),function(ii){class(ii) == "function"}))){
    # If any of the biological rates depend on time  then do simulation with incremental time. If the rates
    # are fixed, the approach is 'waiting time' based, see below
    if (class(spec)!='function'){
      spec_tmp = spec;
      spec <- function(t,n){spec_tmp}}
    if (class(ext) !='function'){
      ext_tmp = ext;
      ext  <- function(t,n){ext_tmp}}
    Times = array(NA,c(n_init*1000,2));
    # Collect times here [start, end]
    Times[1:n_init,1] = 0;
    # Alive = array(0,c(n_init*1000,n_ints));
    # Alive[1:n_init,1]= 1;
    # Obs[1:n_init,1]  = runif(n_init)>samp*1
    # inxtmp = n_init;
    # So we need to integrate forwards in time.
    dxt = min(dt_ints)/(1000*n_init); #temporal resolution
    # Think we collect start and end of all species
    # first
    t_now = 0;
    txmax = n_init; # max total no species for tally
    alive_now = c(1:n_init);
    for (tt in 1:length(dt_ints)){
      while (t_now<sum(dt_ints[1:tt])){
        # For each step, go forwards in time dxt. Calculate probabilities of speciation and extinction in this small window:
        sp_now = max(0,rate2prob((spec(t_now,    length(alive_now)) +
                                    spec(t_now+dxt,length(alive_now)))/2*length(alive_now),dxt));
        ex_now = max(0,rate2prob((ext(t_now,     length(alive_now))  +
                                    ext(t_now+dxt, length(alive_now)))/2*length(alive_now),dxt));
        event <- sample.int(3,1,prob=c(1-(sp_now+ex_now),sp_now,ex_now));
        #if 1 (nothing), 2 speciate, 3 die
        # Random who it happens to
        if (event==1){         # Nothing happened
          t_now = t_now + dxt
        } else if (event==2){  # Speciation event.
          t_now = t_now + dxt;
          alive_now = c(alive_now,txmax+1); # New species
          if (txmax>(-1+dim(Times)[1])){ # If array is full, double it.
            Times = rbind(Times,array(NA,c(dim(Times)[1],dim(Times)[2])))
          }
          Times[txmax+1,1] = runif(1,min=t_now-dxt,max=t_now); # exact time of "birth"
          txmax = txmax+1; #ticker for number of total number of species/entry of next newborn in Times array
        } else if (event==3){ # Extinction event
          t_now = t_now + dxt;
          whichdies <- sample(alive_now,1) # Which species dies
          alive_now = alive_now[-which(alive_now %in% whichdies)]; # remove the species
          Times[whichdies,2] = runif(1,min=t_now-dxt,max=t_now)    # set exact point of death
        }
        if (length(alive_now)==0){ # If all are dead put time to tmax
          t_now = sum(dt_ints);
        }
        if ((sp_now+ex_now)>1/4){ # If probabilities are "high", we are unintentially linearizing the system. Reduce stepsize in time.
          dxt = dxt/2; # half the stepsize
          print(c(sp_now,ex_now,event,t_now/sum(dt_ints),dxt))
        }
        else if (max(sp_now,ex_now)<1e-6){ # If probabilities are low
          dxt = dxt*2; # double the stepsize
          print(c(sp_now,ex_now,event,t_now/sum(dt_ints),dxt)) # This also prints a set of numbers on screen.
        }
      }
    }
    Taxa = Times[1:txmax,]; # scrap remainder of array
  } else {
    # All rates are fixed, waiting time simulations are faster.
    l = spec
    m = ext
    tmax = sum(dt_ints)
    taxa = array(NA,dim=c(max(100,n_init)^2,2))
    taxa[1:n_init,1] = 0;
    alive = 1:n_init;
    ntix = n_init+1; #place next individual here
    # Draging waiting times for all
    # hist(replicate(1e4,which.min(c(min(rexp(n,rate = l)),min(rexp(n,rate = m))))))
    # If death, killing one randomly chosen.
    t = 0;
    while (t<tmax){
      waitingtimes <- cbind(rexp(length(alive),rate = l),rexp(length(alive),rate = m))
      event <- which.min(c(min(waitingtimes[,1]),min(waitingtimes[,2])))
      t = t+min(waitingtimes)
      if (t>=tmax){break}
      # print(event)
      # if (length(alive)==0){
      #   t = tmax
      # }

      if (event==1){
        # birth
        taxa[ntix,1] = t;
        alive = c(alive,ntix); # array of alive species
        ntix = ntix+1;
        if (ntix==dim(taxa)[1]){
          taxa = rbind(taxa,array(NA,dim=c(n_init^2,2)))
        }
      } else if (event==2){
        if (length(alive)>1){
          dies <- sample(alive,1);
        } else {
          dies = alive; # if only 1 lest, sample interprets this as pick from 1:inx
        }
        taxa[dies,2] = t;
        alive <- setdiff(alive,dies);
        if (length(alive)==0){t = tmax}
      }
    }

    Taxa = taxa[1:(ntix-1),]

  }
  # Now we have the durations of all taxa.
  if (!is.null(dim(Taxa))){# if more than one taxa
    Taxa[is.na(Taxa[,2]),2] = sum(dt_ints) # set the ones that are alive to 'end' at end of time
    Taxa[Taxa[,2]==0,2] = sum(dt_ints)     # set the ones that are alive to 'end' at end of time
    Taxa <- Taxa[!Taxa[,1]>(sum(dt_ints)),] # Remove all born after t-max
    Foss <- lapply(1:dim(Taxa)[1],function(ii){sampFosRec(Taxa[ii,1],Taxa[ii,2],samp)});
    FosRec <- array(0,c(sum(sapply(Foss,length)>0),length(dt_ints)));

  } else {
    # If only 1 taxa, still alive?
    if (is.na(Taxa[2])){
      Taxa[2] = sum(dt_ints)}
    Foss <- sampFosRec(Taxa[1],Taxa[2],samp)
    FosRec <- array(0,c(sum(sapply(Foss,length)>0),length(dt_ints)));

  }
  tix = 1;
  for (jj in which(sapply(Foss,length)>0)){
    FosRec[tix,
           rle(sort(sapply(Foss[[jj]],function(ii){which(ii<cumsum(dt_ints))[1]})))$values] <-
      rle(sort(sapply(Foss[[jj]],function(ii){which(ii<cumsum(dt_ints))[1]})))$lengths
    tix = tix+1;
  }
  rownames(FosRec) <- which(sapply(Foss,length)>0)
  out <- list(Taxa = Taxa,Foss = Foss,FosRec=FosRec,dts = dt_ints,
              Spec = spec, Ext = ext, Samp=samp,n_init=n_init)
  attr(out,"class") <- "cmr_simulation";
  return(out)
}


