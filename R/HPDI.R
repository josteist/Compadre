#' Function to extract highest posterior density intervals.
#' HPDI corresponds to the 'shortest' intervals that contain p fraction of the samples.
#' @param x Samples of one parameter or a CMR_fit.
#' @param p probability. Defaults to 0.95.
#' @return returns hpdi and exact probability of the interval, either for the given sample or for each parameter in a CMR_fit.
#' @export
HPDI <- function(x,p = 0.95){
  # Function to extract highest posterior density intervals.
  #
  #
  if (class(x)=="CMR_fit"){
    # This is a fit, then do for 1000? equally spaced samples from last half
    # x = f1[[1]]
    hpdi = array(NA,c(dim(x$Chain)[2],2))
    exact_prob = array(NA,dim(x$Chain)[2])
    for (jj in 1:dim(x$Chain)[2]){
      x1 = x$Chain[-c(1:round(dim(x$Chain)[1]/2)),jj]
      tix <- which.min(sapply(1:(length(x1)-round(length(x1)*p)),
                              function(ii){sort(x1)[round(length(x1)*p)+ii]-sort(x1)[ii]}))
      hpdi[jj,] <- c(sort(x1)[tix],sort(x1)[round(length(x1)*p)+tix])
      exact_prob[jj] <- round(length(x1)*p)/length(x1);

    }
    colnames(hpdi) <- c('lower','upper')
  } else {
    tix <- which.min(sapply(1:(length(x)-round(length(x)*p)),
                            function(ii){sort(x)[round(length(x)*p)+ii]-sort(x)[ii]}))
    hpdi <- c(sort(x)[tix],sort(x)[round(length(x)*p)+tix])
    exact_prob <- round(length(x)*p)/length(x);
  }

  return(list(hpdi = hpdi, prob = exact_prob))
}


