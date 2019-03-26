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

sim_BD_func_v2 <- function(spec=function(t,n){0.1},
                           ext=function(t,n){0.001},
                           samp = function(t){0.3},
                           n_init=100,dt_ints=rep(1,10)){
  # Simulating a birth death process with spec and ext
  # as functions of (time, n_species now). Output
  # is matrix of origination and extinction times. Speciation is budding.

  # Having, spec, ext (and samp) as functions of time(t)
  # (or n(?))
  # fec = 0.17; # prob of species giving rise to new species
  # surv <- function(t) = 1-ext(t); # survival probability.
  # samp = 0.4; # sampling probability
  # n_init = 500; #number of species initially.

  # dt_ints = 50; # list of intervals with duration

  Times = array(NA,c(n_init*1000,2));
  # Collect times here [start, end]
  Times[1:n_init,1] = 0;
  # Alive = array(0,c(n_init*1000,n_ints));
  # Alive[1:n_init,1]= 1;
  # Obs[1:n_init,1]  = runif(n_init)>samp*1
  # inxtmp = n_init;
  # So we need to integrate forwards in time.
  dxt = min(dt_ints)/100*n_init; #temporal resolution
  # Think we collect start and end of all species
  # first
  t_now = 0;
  txmax = n_init; # max total no species for tally
  alive_now = c(1:n_init);
  for (tt in 1:length(dt_ints)){
    print(c(t_now/sum(dt_ints)))
    # looptix = 0; # loopticker
    # nxtix   =  ceiling(sum(dt_ints[1:tt])-t_now )/dxt;
    # chkdxt  <- round(seq(0.1,1,by=9)*nxtix)
    # print(chkdxt)
    while (t_now<sum(dt_ints[1:tt])){
      # looptix = looptix+1;
      # For each interval, go forwards in time.
      # which are alive at tt-1
      # should we do spec(at t) or int(spec(from t0 to t0+dt))
      # we do 'mean', i.e. linearize over dt
      # What is wrong here? It is per lineage!
      # But the prob is not just rate*lineage?
      # We can use the
      sp_now = max(0,rate2prob((spec(t_now,    length(alive_now)) +
                                  spec(t_now+dxt,length(alive_now)))/2*length(alive_now),dxt));
      ex_now = max(0,rate2prob((ext(t_now,     length(alive_now))  +
                                  ext(t_now+dxt, length(alive_now)))/2*length(alive_now),dxt));
      event <- sample.int(3,1,prob=c(1-(sp_now+ex_now),sp_now,ex_now));
      #if 1 (nothing), 2 speciate, 3 die
      # Random who it happens to
      # This is not right, we cannot take the rate and multiply by the number of lines and then transform to probability?
      # We should make a probability of an event for each line, then draw number of lineages in which an event will occur. If this is >1 at the same time, reduce dxt

      # We need to make sure that these event probs are low,
      # if not we are linearizing the system.
      # if (looptix %in% chkdxt){

      #   if (max(sp_now,ex_now)>1/1000){
      #     dxt = dxt/2; # half the stepsize
      #   }
      #   else if (max(sp_now,ex_now)<1e-6){
      #     dxt = dxt*2; # double the stepsize
      #   }
      # # }


      if (event==1){
        # Nothing happened
        t_now = t_now + dxt
      } else if (event==2){
        # Speciation event.
        t_now = t_now + dxt;
        # New species
        alive_now = c(alive_now,txmax+1);
        if (txmax>(-1+dim(Times)[1])){
          # Times = array(NA,c(dim(Times)[1]*2,dim(Times)[2]));
          Times = rbind(Times,array(NA,c(dim(Times)[1],dim(Times)[2])))
        }
        Times[txmax+1,1] = runif(1,min=t_now-dxt,max=t_now);
        txmax = txmax+1;
      } else if (event==3){
        t_now = t_now + dxt;
        whichdies <- sample(alive_now,1)
        alive_now = alive_now[-which(alive_now %in% whichdies)];
        Times[whichdies,2] = runif(1,min=t_now-dxt,max=t_now)#t_now;
      }
      if (length(alive_now)==0){
        t_now = sum(dt_ints);
      }

    if ((sp_now+ex_now)>1/4){
      dxt = dxt/2; # half the stepsize
      print(c(sp_now,ex_now,event,t_now/sum(dt_ints),dxt))


    }
    else if (max(sp_now,ex_now)<1e-6){
      dxt = dxt*2; # double the stepsize
      print(c(sp_now,ex_now,event,t_now/sum(dt_ints),dxt))

    }
    # print(c(sp_now,ex_now,dxt))
    }
    # do<- alive_now;
    #   # which(Alive[,tt-1]==1);
    # speciate <- sample(alive_now,
    #               rbinom(1,length(alive_now),sp_now))
    # # assume speciation occurs before extinction?sp_
    # die      <- sample(alive_now,
    #                    rbinom(1,length(alive_now),ex_now))
  }
  # Now we have the durations of all taxa.
  # return(Taxa = Times[1:txmax,])
  Taxa = Times[1:txmax,];
  Taxa[is.na(Taxa[,2]),2]=sum(dt_ints)
  # Taxa[is.na(Taxa[,2]),2]=t_now; #sum(dt_ints);
  # Or if 0 and array has been enlarged.
  Taxa[Taxa[,2]==0,2]=sum(dt_ints)
  # t_now; #sum(dt_ints);

  Foss <- lapply(1:dim(Taxa)[1],function(ii){sampFosRec(Taxa[ii,1],Taxa[ii,2],
                                                        samp)});
  #
  # which(sapply(Foss,length)>0)  # Actually left a fossil record
  # print(sapply(Foss,length))
  FosRec <- array(0,c(sum(sapply(Foss,length)>0),length(dt_ints)));
  tix = 1;
  for (jj in which(sapply(Foss,length)>0)){
    # sapply(c(0.004,1.2,4.5,4.7,8.1,9),function(ii){which(ii<cumsum(dt_ints))[1]})
    # print(rle(sapply(Foss[[jj]],function(ii){which(ii<cumsum(dt_ints))[1]})))
    # print(Foss[[jj]])
    # FosRec[tix,rle(sapply(Foss[[jj]],function(ii){which(ii<cumsum(dt_ints))[1]}))$values] <-
    # rle(sapply(Foss[[jj]],function(ii){which(ii<cumsum(dt_ints))[1]}))$lengths
    # Changed 180119. The rle doesn't work w/o sorting first
    FosRec[tix,
           rle(sort(sapply(Foss[[jj]],function(ii){which(ii<cumsum(dt_ints))[1]})))$values] <-
      rle(sort(sapply(Foss[[jj]],function(ii){which(ii<cumsum(dt_ints))[1]})))$lengths
    # print(tix/dim(FosRec)[1])

    tix = tix+1;

  }
  # spec=0.1,ext=0,samp = 0.3,
  # n_init=100,dt_ints=rep(1,10)
  rownames(FosRec) <- which(sapply(Foss,length)>0)
  return(list(Taxa = Taxa,Foss = Foss,FosRec=FosRec,dts = dt_ints,
              Spec = spec, Ext = ext, Samp=samp,n_init=n_init))
}

