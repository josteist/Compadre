#' Functino to rescale between time-scales.
#'
#' This function takes the time-points x which belong to scale1 and rescales
#' them to scale 2. Scale1 and scale2 are two sets of ordered boundaries.
#'
#' This function can be used to change a time on Gradstein 2004 to Gradstein 2012 for instance.
#'
#' @param scale1 scale on which the times in x lie
#' @param scale2 scale on which the times x should be changed into
#' @param x times of observations on scale1
#'
#' @return x_new - x on scale 2
#' @export
# Put this function in the package.
rescale_time <- function(scale1,scale2,x){
  # From scale_1 to scale_2. Scales are ordered boundaries.
  # To find new time of point x given in current scale to x
  # in scale_new. Useful for timeseries calibrated on older
  # timescales.
  bin = c(max(which(scale1<=x)),(which(scale1>x))[1])
  bt1 = scale1[bin]
  bt2 = scale2[bin]
  x_new = bt2[2]-(-bt2[1]+bt2[2])*((-x+bt1[2])/(-bt1[1]+bt1[2]))
  # b2-(-t2+b2)*((-x+b1)/(-t1+b1));
  return(x_new)
}
