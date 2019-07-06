#' Function to make simple plot of estimated drivers (and intercepts) in CMR model fit.
#'
#' @param cmrfit A model fit, as output from \link{MCMC_CMR}
#' @export

plotDrivers <- function(cmrfit,nsamp=250,spec=T,ext=T,samp=T){
  # Intercepts are exp(values), rest are as normal. Intercepts are always the first-
  smp <- round(seq(dim(cmrfit$Chain)[1]/2,dim(cmrfit$Chain)[1],length.out=nsamp))
  # one plot for each rate: spec, ext, samp
  # Spec
  olp <- par(no.readonly = TRUE);
  on.exit(par(olp))
  # Genericaly what are optimal no of panels for ANY number. Say max c(2,4)
  # 1 (1,1)
  # 2 (1,2)
  # 3 (2,2)
  # 4 (2,2)
  # 5 (2,3)
  # 6 (2,3)
  # 7-8 (2,4)
  if (spec){
    if (length(cmrfit$Model$inx$specInx)>8){
      par(mfrow=c(2,4))
    } else {
      ptmps <- list(c(1,1),c(1,2),c(2,2),c(2,2),c(2,3),c(2,3),c(2,4),c(2,4));
      par(mfrow=ptmps[[length(cmrfit$Model$inx$specInx)]])
    }

    plot(density(exp(cmrfit$Chain[smp,cmrfit$Model$inx$specInx[1]])),
         main = 'Intercept - mean rate',
         xlab=expression(paste("Speciation rate, ",lambda)),yaxt='n',ylab='',bty='n')
    rug(exp(cmrfit$Chain[smp,cmrfit$Model$inx$specInx[1]]))
    if (length(cmrfit$Model$inx$specInx)>1){
      for (ii in 2:length(cmrfit$Model$inx$specInx)){
        plot(density((cmrfit$Chain[smp,cmrfit$Model$inx$specInx[ii]])),
             main = names(cmrfit$Model$inx$specInx)[ii],
             xlab='Impact on speciation rate',yaxt='n',ylab='',bty='n')
        rug(cmrfit$Chain[smp,cmrfit$Model$inx$specInx[ii]])
        abline(v=0,col='red')

      }
      # Include more info here; ESS? pseudo-p's?
    }
    if (ext){

      # Ext
      if (length(cmrfit$Model$inx$extInx)>8){
        par(mfrow=c(2,4))
      } else {
        ptmps <- list(c(1,1),c(1,2),c(2,2),c(2,2),c(2,3),c(2,3),c(2,4),c(2,4));
        par(mfrow=ptmps[[length(cmrfit$Model$inx$extInx)]])
      }

      plot(density(exp(cmrfit$Chain[smp,cmrfit$Model$inx$extInx[1]])),
           main = 'Intercept - mean rate',
           xlab=expression(paste("Extinction rate, ",mu)),yaxt='n',ylab='',bty='n')
      rug(exp(cmrfit$Chain[smp,cmrfit$Model$inx$extInx[1]]))
      if (length(cmrfit$Model$inx$extInx)>1){
        for (ii in 2:length(cmrfit$Model$inx$extInx)){
          plot(density((cmrfit$Chain[smp,cmrfit$Model$inx$extInx[ii]])),
               main = names(cmrfit$Model$inx$extInx)[ii],
               xlab='Impact on extinction rate',yaxt='n',ylab='',bty='n')
          rug(cmrfit$Chain[smp,cmrfit$Model$inx$extInx[ii]])
          abline(v=0,col='red')

        }
      }
    }
    if (samp){
      # Samp
      if (length(cmrfit$Model$inx$sampInx)>8){
        par(mfrow=c(2,4))
      } else {
        ptmps <- list(c(1,1),c(1,2),c(2,2),c(2,2),c(2,3),c(2,3),c(2,4),c(2,4));
        par(mfrow=ptmps[[length(cmrfit$Model$inx$sampInx)]])
      }

      plot(density(exp(cmrfit$Chain[smp,cmrfit$Model$inx$sampInx[1]])),
           main = 'Intercept - mean rate',
           xlab=expression(paste("Sampling rate, ",rho)),yaxt='n',ylab='',bty='n')
      rug(exp(cmrfit$Chain[smp,cmrfit$Model$inx$sampInx[1]]))
      if (length(cmrfit$Model$inx$sampInx)>1){
        for (ii in 2:length(cmrfit$Model$inx$sampInx)){
          plot(density((cmrfit$Chain[smp,cmrfit$Model$inx$sampInx[ii]])),
               main = names(cmrfit$Model$inx$sampInx)[ii],
               xlab='Impact on sampling rate',yaxt='n',ylab='',bty='n')
          rug(cmrfit$Chain[smp,cmrfit$Model$inx$sampInx[ii]])
          abline(v=0,col='red')

        }
      }
    }

  }
}
