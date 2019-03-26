#' Generate datavectors used in Pradel likelihood calculations.
#'
#' This function takes an observation matrix \emph{Obs} and generates the vectors u,n,v needed to calculate the Pradel likelihood. Used internally by other functions.
#' @param Obs a matrix with size \emph{number of taxa} by \emph{number of intervals}. Each taxa has a row with 0's (unobserved) and 1's (observed) for each interval in the analysis.
#' @return a list containing \emph{Obs}, \emph{u,n,v} and \emph{d}.
#' @export

make_unvd <- function(Obs){
  # Extracting u,n,v,d from an observational matrix.
  # n - number of observed individuals/species over time
  # u - number of invididuals/species seen for the first time at t
  # v - number of individuals/species seen for the last time at
  Obs = 1*Obs>0; # In case somebody inputs other than presence/absences.
  fu = function(x){min(which(x>0))}
  tmpu = apply(Obs,1,fu)
  u = array(0,dim(Obs)[2])
  u[rle(sort(tmpu))$values] = rle(sort(tmpu))$lengths

  n = colSums(Obs)
  vu = function(x){max(which(x>0))}
  tmpv = apply(Obs,1,vu)
  v = array(0,dim(Obs)[2]);
  v[rle(sort(tmpv))$values] = rle(sort(tmpv))$lengths
  d = array(data=0,length(u))
  out = list(Obs,u,n,v,d);
  names(out)<-c('Obs','u','n','v','d');
  return(out)

}
