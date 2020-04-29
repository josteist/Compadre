#' Function to calculate Alroy's three-timer rates of macroevolutionary change
#'
#' @param Obs matrix of observations (taxa by time)
#' @param dts vector of interval duration (defaults to 1)

#' @export
#' @return
Alroy_3timers <- function(Obs,dts=rep(1,dim(Obs)[2])){

  # Making a function to do Alroy's three-timer approach.


  # Need to extract
  # - two-timers: taxa sampled immediately before and within a bin (2t[i])
  #   OR within and immediately after (2t[i+1])
  # (I think this is just poorly written, twotimers are just 2t[i], and we use them diff in calculating lambda/mu below.
  # - three-timers: taxa sapmled immediately before, within and after a
  #   bin (3t[i])
  # - part-timers: taxa samples immediately before and after but not within
  #   the bin (pt[i])

  # Define part-timer sampling prob as
  # P[s] = N3t / (N3t + Npt)
  # ref.
  # https://www.pnas.org/content/105/Supplement_1/11536#sec-12 ALroy PNAS 2008

  # So we want these for each interval (though first and last will
  # not compute)

  twotimers <- c(0,sapply(2:(dim(Obs)[2]-1),
                          function(ii){sum(apply(Obs[,(ii-1):(ii)],1,prod))}),0)
                              # sum(apply(Obs[,(ii):(ii+1)],1,prod))}),0)
  threetimers <- c(0,sapply(2:(dim(Obs)[2]-1),
                            function(ii){sum(apply(Obs[,(ii-1):(ii+1)],1,prod))}),0)
  parttimers <- c(0,sapply(2:(dim(Obs)[2]-1),function(ii){
    sum(Obs[,ii]==0 & apply(Obs[,c(ii-1,ii+1)],1,prod)>0)}),0)

  Ps <- (threetimers)/(threetimers+parttimers)
  Ps2 <- sum(threetimers)/sum(threetimers+parttimers) # global sampling prob
  lambda <- c(NA,sapply(2:dim(Obs)[2],function(ii){
    (log(twotimers[ii+1]/threetimers[ii]) + log(Ps[ii-1]))*(1/dts[ii])}))
  mu <- sapply(1:dim(Obs)[2],function(ii){
    (log(twotimers[ii]/threetimers[ii]) + log(Ps[ii+1]))*(1/dts[ii])})
  names(lambda) <- names(mu)
  lambda_2 <- c(sapply(1:dim(Obs)[2],function(ii){
    (log(twotimers[ii+1]/threetimers[ii]) + log(Ps2))*(1/dts[ii])}))
  mu_2 <- sapply(1:dim(Obs)[2],function(ii){
    (log(twotimers[ii]/threetimers[ii]) + log(Ps2))*(1/dts[ii])})

    out <- list(twotimers = twotimers,threetimers=threetimers,parttimers=parttimers,
              SamplingProb = Ps,mu= mu, lambda=unlist(lambda),
              SamplingProb2 = Ps2,mu2= mu_2, lambda2=unlist(lambda_2))
  return(out)
}
