#' @export
# quick function to check for too probabilities too close to 0 or 1 inside pradel_unvd_gam
mypmin <- function(x,minval){
  if (any(x<minval)){
    x[which(x<minval)] = minval;
  }
  if (any(x>(1-minval))){
    x[which(x>(1-minval))] = 1-minval
  }

  return(x)
}
