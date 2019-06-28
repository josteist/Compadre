#' A function to calculate the probability of extinction for a paraclade.
#' @param m mu - extinction rate
#' @param l lambda - speciation rate
#' @param dt duration of the interval/period
#' @export
#' 
#' 
prob_ext <- function(m,l,dt){(m * (exp((l-m)*dt) - 1) / (l * exp((l-m)*dt)-m))}
