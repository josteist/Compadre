#' Extract and plot macroevolutionary rates from a CMR model fit.
#'
#' @param cmrfit A model fit, as output from \link{MCMC_CMR}
#' @param max_ma (optional) a maximum age for the oldest interval, if not given assumed to be sum of all interval durations in cmrfit
#' @param stages (optional) a data.frame of stages used for the analysis. Values max_ma,min_ma and colors are used to color the plot.
#' @param logax (T/F), use log y-axis for plotting
#' @param draws number of samples drawn from the MCMC chain (single number). Alternatively could be set of values indexing into the MCMC-chain.
#' @param qauntiles gives the upper and lower quantiles for plotting bars around median, and the centre. Defaults to c(0.025,0.5,0.975).
#' @param plot (T/F) if plot figure, defaults to TRUE.
#' @return returns a list of three matrices with sampled rates; $SpecRates, $ExtRates and $SampRates.
#' @export

plotRates <- function(cmrfit,max_ma= sum(cmrfit$Model$dts),stages =NULL,logax = T,draws = 250,quantiles = c(0.025,0.5,0.975),plot=T){
  # have max_ma as input instead of stgs, and possibly a list of colors.
  # thus boundaries are max_ma+c(0,cumsum(cmrfit$Model$dts))
  # also store par(stuff) and return settings to old after plott.
  # plpars <- par();
  olp <- par(no.readonly = TRUE);
  on.exit(par(olp))
  if (!is.null(stages)){
    # If stages are given, take the max as the xmax
    max_ma = max(stages$max_ma)
  }
  bnds <- rev(max_ma+c(0,cumsum(rev(cmrfit$Model$dts))))-sum(cmrfit$Model$dts);
  if (logax){
    tmpl = 'y'
  } else {
    tmpl = ''
  }
  if (length(draws)==1){
    tpl <- ceiling(seq(dim(cmrfit$Chain)[1]/2,dim(cmrfit$Chain)[1],length.out=draws));
  } else {
    tpl = draws
  }
  par(mfrow=c(3,1),mar=c(2,4,.2,1))

  # extracting spc rate quantile
  spec_smps <- exp((sapply(1:length(tpl),function(ii){
    cmrfit$Model$specfun(cmrfit$Chain[tpl[ii],])})))
  tmp <- apply(spec_smps,1,
               mf<-function(k){quantile(k,quantiles)})
  if (is.null(stages)){
    # if stages not given
    plot(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19,
         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',
         log=tmpl,ylab='Speciation rate',
         ylim=c(min(tmp)*0.9,max(tmp)*1.1))
    # axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
  } else {
    # MAKE THE SHADING HERE...
    plot(1000,1000,         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',
         log=tmpl,ylab='Speciation rate',
         ylim=c(min(tmp)*0.9,max(tmp)*1.1))

    for (ii in 1:dim(stages)[1]){
      if (par("ylog") ){
        rect(stages[ii,]$min_ma,10^(par("usr")[3]),
             stages[ii,]$max_ma,10^(par("usr")[4]),
             col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)

      } else {
        rect(stages[ii,]$min_ma,(par("usr")[3]),
             stages[ii,]$max_ma,(par("usr")[4]),
             col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)

      }
    }


    lines(bnds[-c(1,length(bnds))],tmp[2,],
          type="o",lty=0,col='black',pch=19)
    # axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
    #
    # }
  }
#
#   # plot(rev((stgs$max_ma))[-1],tmp[2,],type="o",lty=0,col='black',pch=19,
#   plot(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19,
#        xlim = c(max(bnds)+3,min(bnds-3)),
#        xaxt='n',log=tmpl,ylab='Speciation rate',
#        ylim=c(min(tmp)*0.9,max(tmp)*1.1))
  for (ii in 1:dim(tmp)[2]){
    lines(rep(bnds[ii+1],2),
          tmp[c(1,3),ii])
  }
  abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))

  # extinction quantile
  ext_smps <- exp((sapply(1:length(tpl),function(ii){
    cmrfit$Model$extfun(cmrfit$Chain[tpl[ii],])})))
  tmp <- apply(ext_smps,1,
               mf<-function(k){quantile(k,quantiles)})
  if (is.null(stages)){
    # if stages not given
    plot(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19,
         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',
         log=tmpl,ylab='Extinction rate',
         ylim=c(min(tmp)*0.9,max(tmp)*1.1))
    # axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
  } else {
    # MAKE THE SHADING HERE...
    plot(1000,1000,         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',
         log=tmpl,ylab='Extinction rate',
         ylim=c(min(tmp)*0.9,max(tmp)*1.1))

    for (ii in 1:dim(stages)[1]){
      if (par("ylog") ){
        rect(stages[ii,]$min_ma,10^(par("usr")[3]),
             stages[ii,]$max_ma,10^(par("usr")[4]),
             col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)

      } else {
        rect(stages[ii,]$min_ma,(par("usr")[3]),
             stages[ii,]$max_ma,(par("usr")[4]),
             col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)

      }
    }


    lines(bnds[-c(1,length(bnds))],tmp[2,],
          type="o",lty=0,col='black',pch=19)
    # axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
    #
    # }
  }
  #
  # plot(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19,
  #      xlim = c(max(bnds)+3,min(bnds-3)),
  #      xaxt='n',log=tmpl,ylab='Extinction rate',
  #      ylim=c(min(tmp)*0.9,max(tmp)*1.1))

  for (ii in 1:dim(tmp)[2]){
    # lines(rep(rev((stgs$max_ma))[ii+1],2),
    lines(rep(bnds[ii+1],2),
          tmp[c(1,3),ii])
  }
  abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))

  # sampling quantiles
  smp_smps <- exp((sapply(1:length(tpl),function(ii){
    cmrfit$Model$sampfun(cmrfit$Chain[tpl[ii],])})))
  tmp <- apply(smp_smps,1,
               mf<-function(k){quantile(k,quantiles)})
  if (is.null(stages)){
    # if stages not given
    plot((bnds[-1]+bnds[-length(bnds)])/2,tmp[2,],type="o",lty=0,col='black',pch=19,
         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',
         log=tmpl,ylab='Sampling rate',
         ylim=c(min(tmp)*0.9,max(tmp)*1.1))
    axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
  } else {
    # MAKE THE SHADING HERE...
    plot(1000,1000,         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',
         log=tmpl,ylab='Sampling rate',
         ylim=c(min(tmp)*0.9,max(tmp)*1.1))

         for (ii in 1:dim(stages)[1]){
           if (par("ylog") ){
             rect(stages[ii,]$min_ma,10^(par("usr")[3]),
                  stages[ii,]$max_ma,10^(par("usr")[4]),
                  col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)

           } else {
             rect(stages[ii,]$min_ma,(par("usr")[3]),
                  stages[ii,]$max_ma,(par("usr")[4]),
                  col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)

           }
         }


         lines((bnds[-1]+bnds[-length(bnds)])/2,tmp[2,],
                type="o",lty=0,col='black',pch=19)
         axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
         #
         # }
  }
  for (ii in 1:dim(tmp)[2]){
    lines(rep(((bnds[-1]+bnds[-length(bnds)])/2)[ii],2),
          tmp[c(1,3),ii])
  }
  abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))

  if (!is.null(stages)){
    colbottom(stages)
    # # if bottom colors given. They should be in the order of the intervals.
    # for (ii in 1:length(botcols)){
    #     if (logax){
    #       rect(bnds[ii],10^(par("usr")[3]),bnds[ii+1],10^(par("usr")[3])*1.15,col=toString(botcols[ii]))#stgs[ii,]$color))
    #     } else {
    #       rect(bnds[ii],par("usr")[3],bnds[ii+1],0,col=toString(botcols[ii]))#stgs[ii,]$color))
    #     }
    #
    # }
  }
  #
  # This is for the coloring, try first w/o
  # for (ii in 1:dim(stgs)[1]){
  #   if (logax){
  #     rect(stgs[ii,]$max_ma,10^(par("usr")[3]),stgs[ii,]$min_ma,10^(par("usr")[3])*1.15,col=toString(stgs[ii,]$color))
  #   } else {
  #     rect(stgs[ii,]$max_ma,par("usr")[3],stgs[ii,]$min_ma,0,col=toString(stgs[ii,]$color))
  #   }
  # }
  # abline(v=union(stgs$min_ma,stgs$max_ma),col=rgb(0.1,0.1,0.1,0.1))
  #
  out <- list(SpecRates= spec_smps,ExtRates = ext_smps,SampRates = smp_smps)
  return(invisible(out))
}
