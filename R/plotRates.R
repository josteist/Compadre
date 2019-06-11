#' Plot macroevolutionary rates from a CMR model fit.
#'
#' @param fit_1 A model fit, as output from \link{MCMC_CMR}
#' @param max_ma (optional) a maximum age for the oldest interval.
#' @param botcols (optional) a list of colors for each interval. If supplied the lower part of the lower panel is colored.
#' @export

plotRates <- function(fit_1,max_ma= sum(fit_1$Model$dts),botcols =NULL,logax = T,qs = c(0.025,0.5,0.975)){
  # have max_ma as input instead of stgs, and possibly a list of colors.
  # thus boundaries are max_ma+c(0,cumsum(fit_1$Model$dts))
  # also store par(stuff) and return settings to old after plott.
  # plpars <- par();
  olp <- par(no.readonly = TRUE);
  on.exit(par(olp))
  bnds <- rev(max_ma+c(0,cumsum(rev(fit_1$Model$dts))))-sum(fit_1$Model$dts);
  if (logax){
    tmpl = 'y'
  } else {
    tmpl = ''
  }
  tpl <- ceiling(seq(dim(fit_1$Chain)[1]/2,dim(fit_1$Chain)[1],length.out=250));
  par(mfrow=c(3,1),mar=c(2,4,.2,1))

  # extracting spc rate quantile
  spec_smps <- exp((sapply(1:length(tpl),function(ii){
    fit_1$Model$specfun(fit_1$Chain[tpl[ii],])})))
  tmp <- apply(spec_smps,1,
    mf<-function(k){quantile(k,qs)})

  # plot(rev((stgs$max_ma))[-1],tmp[2,],type="o",lty=0,col='black',pch=19,
  plot(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19,
       xlim = c(max(bnds)+3,min(bnds-3)),
       xaxt='n',log=tmpl,ylab='Speciation rate',
       ylim=c(min(tmp)*0.9,max(tmp)*1.1))
  for (ii in 1:dim(tmp)[2]){
    lines(rep(bnds[ii+1],2),
          tmp[c(1,3),ii])
  }
  abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))

  # extinction quantile
  ext_smps <- exp((sapply(1:length(tpl),function(ii){
    fit_1$Model$extfun(fit_1$Chain[tpl[ii],])})))
  tmp <- apply(ext_smps,1,
    mf<-function(k){quantile(k,qs)})
  plot(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19,
       xlim = c(max(bnds)+3,min(bnds-3)),
       xaxt='n',log=tmpl,ylab='Extinction rate',
       ylim=c(min(tmp)*0.9,max(tmp)*1.1))
  for (ii in 1:dim(tmp)[2]){
    # lines(rep(rev((stgs$max_ma))[ii+1],2),
    lines(rep(bnds[ii+1],2),
          tmp[c(1,3),ii])
  }
  abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))

  # sampling quantiles
  smp_smps <- exp((sapply(1:length(tpl),function(ii){
    fit_1$Model$ratefun[[3]](fit_1$Chain[tpl[ii],])})))
  tmp <- apply(smp_smps,1,
    mf<-function(k){quantile(k,qs)})
  plot((bnds[-1]+bnds[-length(bnds)])/2,tmp[2,],type="o",lty=0,col='black',pch=19,
       xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',
       log=tmpl,ylab='Sampling rate',
       ylim=c(min(tmp)*0.9,max(tmp)*1.1))
  axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
  for (ii in 1:dim(tmp)[2]){
    lines(rep(((bnds[-1]+bnds[-length(bnds)])/2)[ii],2),
    tmp[c(1,3),ii])
  }
  abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))

  if (!is.null(botcols)){
    # if bottom colors given. They should be in the order of the intervals.
    for (ii in 1:length(botcols)){
        if (logax){
          rect(bnds[ii],10^(par("usr")[3]),bnds[ii+1],10^(par("usr")[3])*1.15,col=toString(botcols[ii]))#stgs[ii,]$color))
        } else {
          rect(bnds[ii],par("usr")[3],bnds[ii+1],0,col=toString(botcols[ii]))#stgs[ii,]$color))
        }

    }
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
  return(out)
}
