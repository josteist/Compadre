#' Simple function to add colorbars at the lower end of a plot. Useful for coloring
#' different stages when plotting rates over time.
#' 
#' @param stages - data.frame of stages used. The function uses the min_ma, max_ma and the strings in color. See GSA_timescale.
#' @export

colbottom <- function(stages){
  # This simple function colors the bottom part of the x-axis with
  # colors specified in stages. Structure of stages is as GSA_timescale
  # minimally require, color as string, min_ma and max_ma.
  tmppar <- par()
  par("xpd" = TRUE)
  # xps allows plotting outside the plot region
  # if bottom colors given. They should be in the order of the intervals.
  for (ii in 1:dim(stages)[1]){
    # if (logax){
    if (par("ylog")){
      rect(stages[ii,]$min_ma,10^(par("usr")[3]),
           stages[ii,]$max_ma,10^(par("usr")[3]-1/30*(diff(par("usr")[3:4]))),
           col=toString(stages[ii,]$color))
      #stgs[ii,]$color))
    } else {
      # rect(bnds[ii],par("usr")[3],bnds[ii+1],0,col=toString(botcols[ii]))#stgs[ii,]$color))
      rect(stages[ii,]$min_ma,par("usr")[3],
           stages[ii,]$max_ma,par("usr")[3]-1/30*(diff(par("usr")[3:4])),
           col=toString(stages[ii,]$color))
      
    }
    
    
  }
  suppressWarnings(  par("xpd" = FALSE))
}

