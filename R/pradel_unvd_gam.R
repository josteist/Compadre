#' Internal function used to calculate the likelihood of the Pradel CMR model using the seniority as input.
#'
#' @param gam Seniority in the Pradel model. Either one value or one for each of the n-1 last intervals.
#' @param ext Extinction probabilities. Same size as fec.
#' @param p Sampling probabilities.
#' @param u,n,v,d Values for number observations, each of length equal to number of intervals. u - number of taxa first observed, n - number of taxa observed, v - number of taxa observed for last time. d should be vector of 0's.
#' @param minval Minimum probabilities used for calculations. No probabilities are allowed be be smaller than minval or larger than 1-minval. Defaults to 1e-5.
#' @return  Returns $LogL (log likelihood), and $Ll.numerator and $Ll.denominator.
#' @export


pradel_unvd_gam = function(gam,ext,p,u,n,v,d,minval=1e-5){
  # JOS: 242020: made the checking for values close to 0 and 1 into a funtion()
  # JOS: 120219: changed the indexing so that INPUTTED GAM is length 1 OR t=s-1. Then internally here, the first gam is set to 0.
  #
  # JOS: 180917
  # JOS: 271117 corrected a bug in the denominator. Now identical to Tenan suppmat.
  # JOS: 281117 reparameterized with gam - seniority.
  # Pradel model with data give as
  # u - number of animals/species/taxa obs for the first time at [t]
  # n - number of --- " --- obs at [t]
  # v - number of --- " --- obs for the last time at [t]
  # d - number of --- " --- removed, incase of incomplete sampling. Always 0 for Paleoapplications.
  # These data have length [s], i.e. number of intervals.
  # Thus there are [s-1] transitions.
  # Parameters:
  # fec - fecundity, i.e. speciation probability, length = s-1;
  # ext - extinction probability, length = s-1;
  # p - sampling probability, length = s.
  s = length(u); # number of intervals.
  t = s -1;      # number of transitions.
  if (length(gam)==1){gam = rep(gam,s);}
  if (length(ext)==1){ext = rep(ext,t);}
  if (length(p)==1){p = rep(p,s)};
  # PAd gam,
  gam = c(0,gam);

  # So Speciation probability is not fecundity?
  # Spec = 1- seniority, which is gamma?
  # gamma = phi / (fec + phi)
  # 1 - (1-ext) / (fec + (1-ext))

  # Actually, fec is 's the fecundity rate of the population, so that f(i) is the number of new animals in the population at time i per animal in the population at time i-1, or N(i+1) = N(i) f(i) + N(i) phi(i) '
  # so need not be limited to <1.
  # Changing parameterization.
  # IMplementing some min/max here, so that these do not
  # go to 1 or 0. Say they need to min be 1e-5 and 1-1e-5 resp?
  gam = mypmin(gam,minval);#pmin(1-minval,pmax(gam,minval));
  ext = mypmin(ext,minval);# pmin(1-minval,pmax(ext,minval));
  p   = mypmin(p,minval);# pmin(1-minval,pmax(p,  minval));
  phi = 1 - ext; # Survival is 1 - extinction.
  #phi[t] = 0;    # set to 0, since unestimatable.
  # THis is the confusing with Tenan, they have transition
  # parameters with length (s), which means first or last is
  # not really used. I think that's why they set gamma[1]<- 0 and
  # phi[s] <- 0.
  # Could repar this as 'seniority' - gamma, so spec.prob = 1 - gamma.
  # rho = lambda in Pradel descript.
  # rho = fec + phi; #Growth, fecundity + survival.
  # gam = array(0,s);
  # gam[1] = 0;
  # gam[2:s] = phi[1:t]/rho[1:t];
  rho = phi[1:t]/gam[2:s]
  mu = array(1,s); # Nobody dies because of sampling.
  xi = array(NA,s);
  xi[1] = 1;
  for (ii in 2:s){
    xi[ii] <- (1 - gam[ii]) + (gam[ii]*((1-p[ii-1])))*xi[ii-1];

  }
  chi = array(NA,s);
  chi[s] = 1;
  for (ii in seq(t,1,by=-1)){
    chi[ii] = (1 - phi[ii]) + (phi[ii]*(1 - p[ii+1])*chi[ii+1])
  }
  ll.n <- array(NA,s); #log likelihood numerators.
  # First is slightly different due to summations.
  ll.n[1] <- u[1] * log(xi[1]) +
    (n[1] * log(p[1])) +
    (sum(u[1]) - 0 - n[1]) * (log( 1 - p[1])) +
    (sum(v[2:s]) * log(phi[1])) +
    #(sum(u[2:s]) * log(1 - phi[1]*(1-mu[1]))) +
    (v[1] - d[1])* log(chi[1])
  # ALso last
  ll.n[s] <- u[s]*log(xi[s]) +
    (sum(u[1:t])*(log(phi[t]) - log(rho[t]))) +
    (n[s] * log(p[s])) +
    (sum(u[1:s]) - sum(v[1:t]) - n[s]) * (log(1 - p[s])) +
    (v[s] - d[s]) * log(chi[s])

  # resdt
  for (ii in 2:t){
    ll.n[ii] <- u[ii] * log(xi[ii]) +
      (sum(u[1:(ii-1)])*(log(phi[ii-1]) - log(rho[ii-1]))) +
      n[ii] * log(p[ii]) +
      (sum(u[1:ii]) - sum(v[1:(ii-1)]) - n[ii]) * log( 1 - p[ii]) +
      (sum(v[(ii+1):s]) * log(phi[ii])) +
      #     (n[ii]-d[ii])*log(µ) --> 0
      #     d[ii]*log(1-µ)--> 0
      #      (sum(u[(ii+1):s]) * log(1 - (p[ii]*(1-mu[ii])))) +
      (v[ii]-d[ii]) * log(chi[ii]);
  }
  # Then the denominator. Assuming all(mu) = 1,
  # i.e. sampling does not kill. Makes this a lot easier.
  pm1 = phi; # when µ=1, this collapses to phi
  pm2 = gam; # rho[j-1]/phi[j-1]
  # tmp <- ksi * cumprod(pm1) * cumprod(gam) * p
  ll.d <- sum(u)*log(sum(xi * c(1,cumprod(phi))*  c(sapply(1:(s-2),function(ii){prod(gam[(ii+1):s])}),gam[s],1) * p))
  # Below is old (i.e. before 271117)
  # ll.d <- sum(u) * log(sum(xi * c(1,cumprod(phi))*rev(c(1,cumprod(gam[2:s])))*p));


  LogL <- sum(ll.n) - ll.d;
  out <- list(LogL,ll.n,ll.d,xi,chi,gam,phi,p);
  names(out)<- c('LogL','Ll numerator','Ll','xi','chi','gam','phi','p')
  return(out);


}
