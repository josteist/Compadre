# Fitta fitta fitta fitta - nå fjerna R heile jævla driden.



mod1 <- make_BayesCMR(S1$FosRec>0,dts=S1$dts,RE=c(T,T,T))

x = runif(tmp$npar,min=-.1,max=.1)
tmp$probfun(x)

lfec <- tmp$ratefunc[[1]](x) + tmp$ratefunc[[4]](x)
lext <- tmp$ratefunc[[2]](x) + tmp$ratefunc[[5]](x)
lsmp <- tmp$ratefunc[[3]](x)
dts = S1$dts

mncs <- make_unvd(S1$FosRec>0)

pradel_unvd_gam(    gam = sapply(1:(length(dts)-1),function(ii){ ( exp ( - exp(lext[ii])*(dts[ii])) / (( 1 + exp(lfec[ii]) - exp(lext[ii]))^(dts[ii])))}),
                    ext = sapply(1:(length(dts)-1),function(ii){rate2prob(exp(lext[ii]),dts[ii])}),
                    p   = sapply(1:(length(dts)),function(ii){rate2prob(exp(lsmp[ii]),dts[ii])}),
                    mncs$u,mncs$n,mncs$v,mncs$d)$LogL

tmp$likfun(x)$LogL


# So I guess we should be able to just augment the output from one of these calls...

S1 <- sim_BD_func(spec=function(t,n){0.1 + 0.01*t}, ext = function(t,n){0.052 + 0.009*t},dt_ints= dts)
S2 <- sim_BD_func(spec=function(t,n){0.3}, ext = function(t,n){0.052 + 0.002*t},dt_ints= dts)

# Always call w/o div.dep
mod1 <- make_BayesCMR(S1$FosRec>0,dts=dts,RE=c(T,T,T))
mod2 <- make_BayesCMR(S2$FosRec>0,dts=dts,RE=c(T,T,T))

# Just putting them consex in par array +
mods <- list(mod1,mod2)
nclade = length(mods)
us <- sapply(1:nclade,function(cc){make_unvd(mods[[cc]]$Obs)$u})
ns <- sapply(1:nclade,function(cc){make_unvd(mods[[cc]]$Obs)$n})
vs <- sapply(1:nclade,function(cc){make_unvd(mods[[cc]]$Obs)$v})
ds <- sapply(1:nclade,function(cc){make_unvd(mods[[cc]]$Obs)$d})



# Fcuk we have 'two' interactions, N affecting specrate in selv and others and
# N affecting extrate in self and others.
IntMatSpec = matrix(c(1,0,0,1),nrow=2)
IntMatExt  = matrix(c(1,2,3,4),nrow=2)

mnorm
BigProbf <- function(x){
  tix = 1;
  # Makgin - n_est - functions?
  # making indexes for X-arrays for each clade
  mix <- list();
  for (cc in 1:nclade){
    mix[[cc]] <- c(tix:(tix-1+mods[[cc]]$npar))
    tix = tix+mods[[cc]]$npar
  }
  tmpsp <- (IntMatSpec + tix)*(IntMatSpec>0)
  tmpsp[tmpsp==0] = NA;
  tix <- tix+sum(unique(array(IntMatSpec))>0);
  tmpex <- (IntMatExt  + tix)*(IntMatExt>0)
  tmpex[tmpex==0] = NA;
  # Well, since this only has the first n-1 intervals.
  # I guess that drivers are inside the lfec7lext etc...
  lsmp <- sapply(1:nclade,function(cc){mods[[cc]]$ratefunc[[3]](x[mix[[cc]]])})
  # n_est <- sapply(1:nclade,function(cc){normfun(colSums(mods[[cc]]$Obs>0)/(exp(lsmp[,cc])))})
  # THe above n_est works for ALL intervals, not only the ones - last, for rates.
  # but here we want to take out the last.
  n_est <- sapply(1:nclade,function(cc){mods[[cc]]$n_norm(x[mix[[cc]]])})
  lfec <- sapply(1:nclade,function(cc){mods[[cc]]$ratefunc[[1]](x[mix[[cc]]])})+
    t(x[tmpsp[cc,]]*t(n_est))

  lext <- sapply(1:nclade,function(cc){mods[[cc]]$ratefunc[[2]](x[mix[[cc]]])})+
    t(x[tmpex[cc,]]*t(n_est))
  # THese work if the intervals are the same, if not, then it won't work.
  # Same for the nstuff. Should I continue here OR make it more general.
  # If more general we need to loop
  # Still assuming same dts...

  lls <- sapply(1:nclade,function(cc){
    pradel_unvd_gam(    gam = sapply(1:(length(dts)-1),function(ii){ ( exp ( - exp(lext[ii,cc])*(dts[ii])) / (( 1 + exp(lfec[ii,cc]) - exp(lext[ii,cc]))^(dts[ii])))}),
                        ext = sapply(1:(length(dts)-1),function(ii){rate2prob(exp(lext[ii,cc]),dts[ii])}),
                        p   = sapply(1:(length(dts)),function(ii){rate2prob(exp(lsmp[ii,cc]),dts[ii])}),
                        us[,cc],ns[,cc],vs[,cc],ds[,cc])$LogL})
  # PRob of random effects for each para in ech chain. Perhasp make a 'xinx min and max' for each model?

  erps <- sapply(1:nclade,function(cc){
    sapply(1:length(erix[[cc]]),function(jj){
      dnorm(mods[[cc]]$reix[[jj]])
    }



    # BUT now we need all the 'n_ests', since we can possibly have different lengths
    # of the intervals.
    n_estf <- function(x){mncs$n/(rate2prob(exp(tmp$ratefunc[[3]](x)),dts))}
    n_estf(x)

    # So we can loop over clades and call make_BayesCMR and use
    # ratefunctions [1:5] of x
    #
    # need reindexig of the X

    # so external full probability for S1;
    prob1 <- function(x,mncs,mod1){
      lf <- mod1$ratefunc[[1]](x) + mod1$ratefunc[[4]](x)

    }
