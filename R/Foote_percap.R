#' Function to calculate Foote's rates of macroevolutionary change.
#'
#' @param Obs matrix of observations (species by time)
#' @param dts vector of interval durations
#'
#' @export
#' @return p_hat,q_hat, N_b_t, N_t,N_b
# Foote rate's function?
Foote_percap <- function(Obs,dts){
  ii = 2; #interval
  tmpObs = array(0,dim(Obs));
  for (tt in 1:dim(Obs)[1]){
    tmpObs[tt,min(which(Obs[tt,]==1)):max(which(Obs[tt,]==1))] = 1
  }

  # TmpObs has gaps filled
  # N_t crosses top boundary, last is 0
  N_t <- c(colSums(sapply(1:(dim(tmpObs)[2]-1),function(ii){
    sapply(1:dim(tmpObs)[1],function(tt){
      # Crosses top boundary and is here?
      # Is observed anytime after ii and at or before ii
      (sum(tmpObs[tt,(ii+1):dim(tmpObs)[2]])>0)*(tmpObs[tt,ii])
      # is found after
      #is foudn here OR before
    })
  })),0)


  # Crosses bottom boundary
  # First is 0, calc for intervals 2 and to ii
  N_b <- c(0,colSums(sapply(2:(dim(tmpObs)[2]),function(ii){
    sapply(1:dim(tmpObs)[1],function(tt){
      # Crosses bottom boundary and is here?
      #
      (sum(tmpObs[tt,(1:(ii-1))])>0)*(sum(tmpObs[tt,ii]))
      # is found after
      #is foudn here OR before
    })
  })))


  # N_bt crosses both bounadaries
  # first and last is 0, calc from 2 to second last int.
  N_b_t <- c(0,colSums(sapply(2:(dim(tmpObs)[2]-1),function(ii){
    sapply(1:dim(tmpObs)[1],function(tt){
      # Crosses bottom boundary and top boundary
      (sum(tmpObs[tt,1:(ii-1)])>0)*((sum(tmpObs[tt,(ii+1):dim(tmpObs)[2]])>0))
      # is found after
      #is foudn here OR before
    })
  })),0)

  # p_hat <-   -log(N_b_t/N_t) / deltaT
  # q_hat <-   -log(N_b_t/N_b) / deltaT
  p_hat <- -log(N_b_t/N_t)/dts
  q_hat <- -log(N_b_t/N_b)/dts
       out = list(p_hat = p_hat,
                  q_hat = q_hat,
                  N_bt  = N_b_t,
                  N_t   = N_t,
                  N_b   = N_b)
       return(out)
}
