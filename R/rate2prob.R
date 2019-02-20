#' A simple function to translate a rate into a probability.
#' @param x rate
#' @param dt interval
#' @return probability.
#' @export
rate2prob <- function(x,dt){1-exp(-x*dt)}
