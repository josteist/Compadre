#' @export
normfun <- function(x){
  # Simple normalizing function for any series.
  return((x-mean(x))/sd(x))
}
