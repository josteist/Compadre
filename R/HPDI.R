#' Function to extract highest posterior density intervals.
#' HPDI corresponds to the 'shortest' intervals that contain p fraction of the samples.
#' @param x Samples of one parameter.
#' @param p probability. Defaults to 0.95.
#' @return returns hpdi and exact probability of the interval
#' @export
HPDI <- function(x,p = 0.95){
  # Function to extract highest posterior density intervals.
  #
  #
  tix <- which.min(sapply(1:(length(x)-round(length(x)*p)),
                          function(ii){sort(x)[round(length(x)*p)+ii]-sort(x)[ii]}))
  hpdi <- c(sort(x)[tix],sort(x)[round(length(x)*p)+tix])
  exact_prob <- round(length(x)*p)/length(x);


  return(list(hpdi = hpdi, prob = exact_prob))
}


