#' A quick function to plot chains and calc ESS's for drivers and mean effects.
#'
#' @param cmrfit A fitted model using MCMC_CMR
#' @export


checkchains <- function(cmrfit){
  olp <- par(no.readonly = TRUE);
  on.exit(par(olp))

  if (is.null(cmrfit$Model$Clade1Mod)){
    # if single clade model
    # Four panels:
    # histogram of ESS for ALL rates (possibly with Re's in diff colors?)
    # then matplotchains for mean/drivers with ESS's written in plot?
    esstmp <- ESS(cmrfit)
    brks <- seq(floor(min(esstmp/10)-1)*10,ceiling(max(esstmp)/10+1)*10,length.out=21)
    par(mfrow=c(2,2))
    # cls = c('black','red','green','lightblue')
    cls = c(rgb(0,0,0,.5),   rgb( 1,0,0,.5),   rgb( 0,1,0,.5))
    matplot(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),1:3],type="l",col = cls,lty=1,xlab='iteration',ylab='means',main='')
    if (length(unlist(cmrfit$Model$inx[1:3]))>3){
      # if any drivers
      matplot(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),setdiff(c(unlist(cmrfit$Model$inx[1:3])),1:3)],type="l",lty=1,
              col = c(rep(cls[1],length(cmrfit$Model$inx$specInx)-1),
                      rep(cls[2],  length(cmrfit$Model$inx$extInx)-1),
                      rep(cls[3],  length(cmrfit$Model$inx$sampInx)-1)),ylab='drivers',xlab='iteration')
    } else {
      plot.new()
    }
    cls = c(rgb(200/255,  20/255, 20/255,1),
            rgb( 20/255, 200/255, 20/255,1),
            rgb( 20/255,  20/255,200/255,1))
    hist(esstmp,breaks=brks,col=cls[3],
         main='Effective sample sizes',xlab='ESS',ylab='# parameters')
    hist(esstmp[c(unlist(cmrfit$Model$inx[1:3]))],breaks=brks,col=cls[2],add = TRUE)
    hist(esstmp[1:3],breaks=brks,col=cls[1],add = TRUE)
    if (max(unlist(cmrfit$Model$inx))){
    plot.new()
    legend(x=c(0,1),y=c(0,1),c('Means','Drivers','REs'),fill=cls,bty='n',cex=0.8)
    }
  } else {
    # two clade model
    # 8 panels?
    # Well, perhaps 6: histogram (all, color means, color drivers, color RE's for clades indep
    # then matplot of pars[1:3], then of all drivers for each model (color for which driver))
    # return(ESSs=coda::effectiveSize(coda::as.mcmc(f1$Chain[-c(1:dim(f1$Chain)[2]/2),])))
    esstmp <- ESS(cmrfit)
    brks <- seq(floor(min(esstmp/10)-1)*10,ceiling(max(esstmp)/10+1)*10,length.out=21)
    par(mfrow=c(2,3))
    # cls = c('black','red','green','lightblue')
    cls = c(rgb(0,0,0,.5),   rgb( 1,0,0,.5),   rgb( 0,1,0,.5))
    matplot(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),1:3],type="l",col = cls,lty=1,xlab='iteration',ylab='means',main='clade 1')
    if (length(unlist(cmrfit$Model$inx$inx1[1:3]))>3){
      # if any drivers
      matplot(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),setdiff(c(unlist(cmrfit$Model$inx$inx1[1:3])),1:3)],type="l",lty=1,
              col = c(rep(cls[1],length(cmrfit$Model$inx$inx1$specInx)-1),
                      rep(cls[2],  length(cmrfit$Model$inx$inx1$extInx)-1),
                      rep(cls[3],  length(cmrfit$Model$inx$inx1$sampInx)-1)),ylab='drivers',xlab='iteration')
    } else {
      plot.new()
    }
    cls = c(rgb(200/255,  20/255, 20/255,1),
            rgb( 20/255, 200/255, 20/255,1),
            rgb( 20/255,  20/255,200/255,1))
    hist(esstmp[c(unlist(cmrfit$Model$inx$inx1))],breaks=brks,col=cls[3],
         main='Effective sample sizes',xlab='ESS',ylab='# parameters')
    hist(esstmp[c(unlist(cmrfit$Model$inx$inx1[1:3]))],breaks=brks,col=cls[2],add = TRUE)
    hist(esstmp[1:3],breaks=brks,col=cls[1],add = TRUE)
    legend("topright",c('Means','Drivers','REs'),fill=cls)



    # Clade 2
    cls = c(rgb(0,0,0,.5),   rgb( 1,0,0,.5),   rgb( 0,1,0,.5))

    matplot(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),min(unlist(cmrfit$Model$inx$inx2))+c(0:2)],type="l",col = cls,lty=1,xlab='iteration',ylab='means',main='clade 2')
    if (length(unlist(cmrfit$Model$inx$inx2[1:3]))>3){
      # if any drivers

      matplot(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),setdiff(c(unlist(cmrfit$Model$inx$inx2[1:3])),min(unlist(cmrfit$Model$inx$inx2))+c(0:2))],type="l",lty=1,
              col = c(rep(cls[1],length(cmrfit$Model$inx$inx2$specInx)-1),
                      rep(cls[2],  length(cmrfit$Model$inx$inx2$extInx)-1),
                      rep(cls[3],  length(cmrfit$Model$inx$inx2$sampInx)-1)),ylab='drivers',xlab='iteration')
    } else {
      plot.new();
    }
    cls = c(rgb(200/255,  20/255, 20/255,1),
            rgb( 20/255, 200/255, 20/255,1),
            rgb( 20/255,  20/255,200/255,1))
    hist(esstmp[c(unlist(cmrfit$Model$inx$inx2))],breaks=brks,col=cls[3],
         main='Effective sample sizes',xlab='ESS',ylab='# parameters')
    hist(esstmp[c(unlist(cmrfit$Model$inx$inx2[1:3]))],breaks=brks,col=cls[2],add = TRUE)
    hist(esstmp[min(unlist(cmrfit$Model$inx$inx2))+c(0:2)],breaks=brks,col=cls[1],add = TRUE)
    legend("topright",c('Means','Drivers','REs'),fill=cls)

  }
}
