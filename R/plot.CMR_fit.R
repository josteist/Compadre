#' @export
plot.CMR_fit <- function(fit,drawplot=TRUE,...){
  if (drawplot){
    olp <- par(no.readonly = TRUE);
    on.exit(par(olp))
  }
  # simple wrapper to plot both the rates and the drivers.
  # a switch, if only three rates (i.e. no RE's or drivers, then just plot the density of the posterior)
  if (is.null(fit$Model$Clade1Mod)){ #If not interacting clades
    if (dim(fit$Chain)[2]==3){
      # only three rates
      if (fit$Model$modeltype %in% c("I","II","III","IV")){
        tmp <- exp(fit$Chain[-c(1:dim(fit$Chain)[1]/2),]);
        if (drawplot){
          d1 <- density(exp(fit$Chain[-c(1:dim(fit$Chain)[1]/2),1]),from=0)
          d2 <- density(exp(fit$Chain[-c(1:dim(fit$Chain)[1]/2),2]),from=0)
          d3 <- density(exp(fit$Chain[-c(1:dim(fit$Chain)[1]/2),3]),from=0)
          plot(d1,col=rgb(0.1,0.1,0.8,0.6),
               xlim=c(1e-5,ceiling(10*exp(max(fit$Chain[-c(1:dim(fit$Chain)[1]/2)])))/10),main='Macroevolutionary rates',ylab='',
               xlab='Rate',yaxt='n',bty='n',
               ylim=c(0,max(c(d1$y,d2$y,d3$y))))
          polygon(c(1e-8,d1$x,1e8),
                  c(1e-8,d1$y,1e-8),col=rgb(0.1,0.1,0.8,0.3),border=rgb(0.1,0.1,0.8,0.6))
          polygon(c(1e-8,d2$x,1e8),
                  c(1e-8,d2$y,1e-8),col=rgb(0.8,0.1,0.1,0.3),border=rgb(0.8,0.1,0.1,0.6))
          polygon(c(1e-8,d3$x,1e8),
                  c(1e-8,d3$y,1e-8),col=rgb(0.1,0.6,0.1,0.3),border=rgb(0.1,0.6,0.1,0.6))

          legend("topright",legend=c("Speciation rate","Extinction rate","Sampling rate"),
                 fill= c(rgb(0.1,0.1,0.8,0.3),
                         rgb(0.8,0.1,0.1,0.3),
                         rgb(0.1,0.8,0.1,0.3)),
                 border= c(rgb(0.1,0.1,0.8,0.6),
                           rgb(0.8,0.1,0.1,0.6),
                           rgb(0.1,0.8,0.1,0.6)))
        }
        out <- list(SpecRates = tmp[,1],ExtRates=tmp[,2],SampRates = tmp[,3]);
        } else {
          tmp <- 1/(1+exp(-fit$Chain[-c(1:dim(fit$Chain)[1]/2),]))
          if (drawplot){
            d1 <- density(1/(1+exp(-(fit$Chain[-c(1:dim(fit$Chain)[1]/2),1]))),from=0,to=1)
            d2 <- density(1/(1+exp(-(fit$Chain[-c(1:dim(fit$Chain)[1]/2),2]))),from=0,to=1)
            d3 <- density(1/(1+exp(-(fit$Chain[-c(1:dim(fit$Chain)[1]/2),3]))),from=0,to=1)
            plot(d1,col=rgb(0.1,0.1,0.8,0.6),
                 xlim=c(0,1),main='Macroevolutionary probabilities',ylab='',
                 xlab='Probability',yaxt='n',bty='n',
                 ylim=c(0,max(c(d1$y,d2$y,d3$y))))
            polygon(c(1e-8,d1$x,1e8),
                    c(1e-8,d1$y,1e-8),col=rgb(0.1,0.1,0.8,0.3),border=rgb(0.1,0.1,0.8,0.6))
            polygon(c(1e-8,d2$x,1e8),
                    c(1e-8,d2$y,1e-8),col=rgb(0.8,0.1,0.1,0.3),border=rgb(0.8,0.1,0.1,0.6))
            polygon(c(1e-8,d3$x,1e8),
                    c(1e-8,d3$y,1e-8),col=rgb(0.1,0.6,0.1,0.3),border=rgb(0.1,0.6,0.1,0.6))

            legend("topright",legend=c("Seniority probability","Extinction probabiliry","Sampling probability"),
                   fill= c(rgb(0.1,0.1,0.8,0.3),
                           rgb(0.8,0.1,0.1,0.3),
                           rgb(0.1,0.8,0.1,0.3)),
                   border= c(rgb(0.1,0.1,0.8,0.6),
                             rgb(0.8,0.1,0.1,0.6),
                             rgb(0.1,0.8,0.1,0.6)))
          }
          out <- list(SeniorityProbs = tmp[,1],
                      ExtinctionProbs  = tmp[,2],
                      SamplingProbs = tmp[,3])
        }
        tmp <- out;
      } else {
        tmp <- invisible(plotRates(fit,drawplot=drawplot,...))
      }
      # Also input a plot of the drivers here.
    } else {
      # Interacting clades.
      cat('Plotting for interacting clades is not developed yet. Figure it out for yourself!')
      tmp <- invisible(plotRates(fit,drawplot=drawplot,...))
    }
    return(invisible(tmp))
  }
