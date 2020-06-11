#' Plots simulation output in Compadre.
#'
#' Produces simple figure with the true number of taxa and the unique number of species observed within each interval.
#' @param Sim1 Output from sim_bd_func
#' @export
plot.cmr_simulation <- function(Sim1,drawplot=T){
  x <- c(0,sort(Sim1$Taxa[Sim1$Taxa>0 & Sim1$Taxa<max(Sim1$Taxa)]))
  y <- c(sum(Sim1$Taxa[,1]==0),
         sum(Sim1$Taxa[,1]==0) + cumsum(sapply(2:length(x),function(ii){(x[ii] %in% Sim1$Taxa[,1])*1 - (x[ii] %in% Sim1$Taxa[,2])})))
  if (drawplot){
  plot(x,y,type="l",ylim=c(0,max(y)*1.1),xlab='Time',ylab='#Species / Fossils')
  abline(v = c(0,cumsum(Sim1$dts)),col='grey')
  lines((c(0,cumsum(Sim1$dts[-length(Sim1$dts)])) + cumsum(Sim1$dts))/2,colSums(Sim1$FosRec>0),type="o")
  }
  return(list(time=x,no_spec=y))
}
