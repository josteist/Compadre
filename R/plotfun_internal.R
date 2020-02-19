#' internal plottingfunction
#' export
plotfun_internal <- function(spec_smps,ext_smps,smp_smps,bnds,tmpl,stages=NULL,drawplot=TRUE,quantiles=c(0.025,0.5,0.975),main=''){
  # extracting spc rate quantile
  tmp <- apply(spec_smps,1, mf<-function(k){quantile(k,quantiles)})
  if (drawplot){
    if (is.null(stages)){
      # if stages not given
      plot(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19,xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',log=tmpl,ylab='Speciation rate',ylim=c(min(tmp)*0.9,max(tmp)*1.1),xlab='',main=main)
    } else {
      # MAKE THE SHADING HERE...
      plot(1000,1000,         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',       log=tmpl,ylab='Speciation rate',       ylim=c(min(tmp)*0.9,max(tmp)*1.1),xlab='',main=main)

      for (ii in 1:dim(stages)[1]){
        if (par("ylog") ){
          rect(stages[ii,]$min_ma,10^(par("usr")[3]),           stages[ii,]$max_ma,10^(par("usr")[4]),           col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)
        } else {
          rect(stages[ii,]$min_ma,(par("usr")[3]),           stages[ii,]$max_ma,(par("usr")[4]),           col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)
        }
      }
    }
    lines(bnds[-c(1,length(bnds))],tmp[2,], type="o",lty=0,col='black',pch=19)
    for (ii in 1:dim(tmp)[2]){
      lines(rep(bnds[ii+1],2),tmp[c(1,3),ii])
    }
    abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))

  }

  # extinction quantile
  tmp <- apply(ext_smps,1, mf<-function(k){quantile(k,quantiles)})
  if (drawplot){
    if (is.null(stages)){
      # if stages not given
      plot(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19,xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',log=tmpl,ylab='Extinction rate',ylim=c(min(tmp)*0.9,max(tmp)*1.1),xlab='')
    } else {
      # MAKE THE SHADING HERE...
      plot(1000,1000,         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',log=tmpl,ylab='Extinction rate',ylim=c(min(tmp)*0.9,max(tmp)*1.1),xlab='')

      for (ii in 1:dim(stages)[1]){
        if (par("ylog") ){
          rect(stages[ii,]$min_ma,10^(par("usr")[3]),             stages[ii,]$max_ma,10^(par("usr")[4]),col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)
        } else {
          rect(stages[ii,]$min_ma,(par("usr")[3]),stages[ii,]$max_ma,(par("usr")[4]),col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)
        }
      }
      lines(bnds[-c(1,length(bnds))],tmp[2,],type="o",lty=0,col='black',pch=19)

    }

    for (ii in 1:dim(tmp)[2]){
      lines(rep(bnds[ii+1],2),tmp[c(1,3),ii])
    }
    abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))
  }

  # sampling quantiles
  tmp <- apply(smp_smps,1,mf<-function(k){quantile(k,quantiles)})
  if (drawplot){
    if (is.null(stages)){
      # if stages not given
      plot((bnds[-1]+bnds[-length(bnds)])/2,tmp[2,],type="o",lty=0,col='black',pch=19,xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',log=tmpl,ylab='Sampling rate',ylim=c(min(tmp)*0.9,max(tmp)*1.1),xlab='')
      axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
    } else {
      # MAKE THE SHADING HERE...
      plot(1000,1000,         xlim = c(max(bnds)+3,min(bnds-3)),xaxt='n',log=tmpl,ylab='Sampling rate',ylim=c(min(tmp)*0.9,max(tmp)*1.1),xlab='')

      for (ii in 1:dim(stages)[1]){
        if (par("ylog") ){
          rect(stages[ii,]$min_ma,10^(par("usr")[3]),stages[ii,]$max_ma,10^(par("usr")[4]),col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)
        } else {
          rect(stages[ii,]$min_ma,(par("usr")[3]),stages[ii,]$max_ma,(par("usr")[4]),col=toString(paste(stages[ii,]$color,'44',sep='')),border=NA)
        }
      }
      lines((bnds[-1]+bnds[-length(bnds)])/2,tmp[2,],type="o",lty=0,col='black',pch=19)
      axis(1,at=bnds[round(seq(1,length(bnds),length.out=5))])
    }
    for (ii in 1:dim(tmp)[2]){
      lines(rep(((bnds[-1]+bnds[-length(bnds)])/2)[ii],2),tmp[c(1,3),ii])
    }
    abline(v=bnds,col=rgb(0.1,0.1,0.1,0.1))
  }
}
