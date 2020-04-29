#' Function to calculate Foote's rates of macroevolutionary change.
#'
#' @param Obs matrix of observations (taxa by time)
#' @param dts vector of interval durations (defaults to 1)
#'
#' @export
#' @return lambda,mu, bottomboundcros, topboundcros, rangethrough
# Foote rate's function?

Foote_percap <- function(Obs,dts=rep(1,dim(Obs)[2])){
  firstlast <- sapply(1:dim(Obs)[1],function(ii){range(which(Obs[ii,]>0))})
  bottomboundcros <- sapply(1:dim(Obs)[2],function(ii){
    sum(firstlast[2,]==ii & firstlast[1,]<ii)})
  topboundcros    <- sapply(1:dim(Obs)[2],function(ii){
    sum(firstlast[1,]==ii & firstlast[2,]>ii)})
  rangethrough <- sapply(1:dim(Obs)[2],function(ii){
    sum(firstlast[1,]<ii & firstlast[2,]>ii)})

  lambda <- -log(rangethrough/(rangethrough+topboundcros))*(1/dts)
  mu <- -log(rangethrough/(rangethrough+bottomboundcros))*(1/dts)
  out = list(lambda  = lambda,
             mu = mu,
             bottomboundcros=bottomboundcros ,
             topboundcros = topboundcros,
             rangethrough = rangethrough)
}
