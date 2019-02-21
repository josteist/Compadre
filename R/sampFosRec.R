#' Function to sample fossil from a lineage.
#' 
#' A simple function that yields a 'fossil record' of a lineage. Frequently used in order to simulate
#' fossil records for testing models, perform adequacy analyses or just for fun.
#' 
#' @param orig time of origination of the lineage
#' @param ext time of extinction of a lineage
#' @param samp how fossils are sampled. If given as a number it is interpreted as a proper rate (with units fossils / unit time). Can also be given as a function of time taking argument (t). Note that t is in absolute terms and not from time of origin of the lineage.
#' @return an array of times at which the lineage is sampled
#' @examples 
#' Fos <- sampFosRec(orig = 2,ext = 12,samp = function(t){0.1+0.01*t})
#' Here a lineage originated at time 2 and went extinct at time 12 is sampled according to a linearly increasing samplin rate.
#' @export


# Samp (if function) has (absolute time, age of lineage) as
# input for sampling rate function
sampFosRec <- function(orig,ext,samp){
  if (is.function(samp)){
    # sample as function.
    # Samp only has t OR duration as input?
    # t for absolute time (related to orig/ext)
    # dt for time since origination (distance to orig)
    t_now = orig;
    dxt   = (-orig+ext)/1000; # 1000 bins for each longevity?
    Fos = c();
    while (t_now<ext){
      if (runif(1)<max(0,samp(t_now,t_now-orig))*dxt){
        Fos = c(Fos,t_now)}
      t_now <- t_now + dxt;
    }
  } else {
    # Sampling is fixed for entire duration of lineage
    # Number is poisson distributed, and times are unif distr
    Fos <- runif(rpois(1,samp*(ext-orig)),min=orig,max=ext);
  }
  return(Fos)
}