#' Function transform a probability into a rate for an interval
#' @param p probability of an event
#' @param dt duration of interval
#' @return The function returns the transformed rate -log(1-p)/dt
#' @export
prob2rate <- function(p,dt){
  -log(1-p)/dt
}