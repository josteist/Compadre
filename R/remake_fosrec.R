#' Function to resample a fossil record from a simulated clade.
#' @param sim sets the type parameterization of the rates of speciation and extinction. 'I' and 'II' specifies that the probabilities of extinction and seniority are calculated based on actual durations of intervals, whereas 'III' and 'IV' uses the differences between mid-points of intervals. 'I' and 'III' uses a basic implementation of probability of extinction as 1 - exp(-rate * dt). Models 'II' and 'IV' use a more detailed functional relationship between rate and probability, following Raup (1984):  (µ*((exp((lambda-µ)*dt))-1))/(lambda*((exp((lambda-µ)*dt)))-µ). Modeltype 'V' implements a temporally uninformed model, estimating the probabilities directly (or rather their logit).
#' @param samp function detailing sampling rate with arguments (t,n)
#' @return The function returns an object of class CMR_model. This most important output is out$probfun which is a function for the posterior. This output can be fed directly into the sampler \link{MCMC_CMR}. It also returns \emph{Obs} and \emph{dts} as well as the full call (exlcluding default settings), \emph{call}
#' @export

# Making a function to 'resample' a full set of taxa from a simulation
remake_fosrec <- function(sim = sim, samp = function(t,n){0.01},newdts=NULL){
  sim_out = sim;
  Foss <- lapply(1:dim(sim$Taxa)[1], function(ii) {
    sampFosRec(sim$Taxa[ii, 1], sim$Taxa[ii, 2], samp)
  })
  if (is.null(newdts)){
    FosRec <- array(0, c(sum(sapply(Foss, length) > 0), length(sim$dts)))
    tix = 1
    for (jj in which(sapply(Foss, length) > 0)) {
      FosRec[tix, rle(sort(sapply(Foss[[jj]], function(ii) {
        which(ii < cumsum(sim$dts))[1]
      })))$values] <- rle(sort(sapply(Foss[[jj]], function(ii) {
        which(ii < cumsum(sim$dts))[1]
      })))$lengths
      tix = tix + 1
    }
  } else {
    # new temporal divisions; assuming it covers the same period
    FosRec <- array(0, c(sum(sapply(Foss, length) > 0), length(sim$dts)))
    tix = 1
    for (jj in which(sapply(Foss, length) > 0)) {
      FosRec[tix, rle(sort(sapply(Foss[[jj]], function(ii) {
        which(ii < cumsum(newdts))[1]
      })))$values] <- rle(sort(sapply(Foss[[jj]], function(ii) {
        which(ii < cumsum(newdts))[1]
      })))$lengths
      tix = tix + 1
    }
    sim_out$dts = newdts
  }
  sim_out$Foss = Foss;
  sim_out$FosRec = FosRec;
  sim_out$Samp = samp;

  return(sim_out)
}
