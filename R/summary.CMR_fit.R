#' @export
# make a 'summary' function for a fit.

summary.CMR_fit <- function(cmrfit,nsamp = 1e4){
  # quick stats function for a named variable.
  # Ideally this should also store the output as a matrix or dframe
  # Output mean,p2.5,p50,p97.5,p>/<(cmrfit$Chain)[1]/2,
  if (is.null(cmrfit$Model$Clade1Mod)){
    # Single clade model
    smp <- round(seq(dim(cmrfit$Chain)[1]/2,
                     dim(cmrfit$Chain)[1],length.out=nsamp))
    
    tmpspec <- array(NA,c(length(cmrfit$Model$inx$specInx),6))
    tmpnam <- c()
    for (ii in 1:length(cmrfit$Model$inx$specInx)){
      tmp <- cmrfit$Chain[smp,cmrfit$Model$inx$specInx[ii]];
      tmpnam[[ii]] <- noquote(names(cmrfit$Model$inx$specInx)[[ii]]);
      if (ii==1){
        tmp = exp(tmp);
        tmpnam[[ii]] = 'Overall rate';
      }
      tmpspec[ii,] = c(mean(tmp),           quantile(tmp,c(0.025,0.5,0.975)),sum(tmp*(sign(mean(tmp)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$specInx[ii]]))
    }
    colnames(tmpspec) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmpspec) <- tmpnam;
    ### Extinction
    tmpext <- array(NA,c(length(cmrfit$Model$inx$extInx),6))
    tmpnam <- c()
    for (ii in 1:length(cmrfit$Model$inx$extInx)){
      tmp <- cmrfit$Chain[smp,cmrfit$Model$inx$extInx[ii]];
      tmpnam[[ii]] <- noquote(names(cmrfit$Model$inx$extInx)[[ii]]);
      if (ii==1){
        tmp = exp(tmp);
        tmpnam[[ii]] = 'Overall rate';
      }
      tmpext[ii,] = c(mean(tmp),           quantile(tmp,c(0.025,0.5,0.975)),sum(tmp*(sign(mean(tmp)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$extInx[ii]]))
    }
    colnames(tmpext) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmpext) <- tmpnam;
  ## == SAMPLING
    tmpsamp <- array(NA,c(length(cmrfit$Model$inx$sampInx),6))
    tmpnam <- c()
    for (ii in 1:length(cmrfit$Model$inx$sampInx)){
      tmp <- cmrfit$Chain[smp,cmrfit$Model$inx$sampInx[ii]];
      tmpnam[[ii]] <- noquote(names(cmrfit$Model$inx$sampInx)[[ii]]);
      if (ii==1){
        tmp = exp(tmp);
        tmpnam[[ii]] = 'Overall rate';
      }
      tmpsamp[ii,] = c(mean(tmp),           quantile(tmp,c(0.025,0.5,0.975)),sum(tmp*(sign(mean(tmp)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$sampInx[ii]]))
    }
    colnames(tmpsamp) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmpsamp) <- tmpnam;
    writeLines("\t \t \t=== Speciation rate parameters ===")
    print(tmpspec,digits=3)
    writeLines("\t \t \t=== Extinction rate parameters ===")
    print(tmpext,digits=3)
    writeLines("\t \t \t=== Sampling rate parameters ===")
    print(tmpsamp,digits=3)
    out = list(SpecStat=tmpspec,
               ExtStat =tmpext,
               SampStat=tmpsamp);
    # rownames(tmpspec) <- c('Overall rate',noquote(names(cmrfit$Model$inx$specInx)[[2]]))
    # format(tmpspec,digits=3)
    return(invisible(out));
  } else {
    # Two clade model...
    # Single clade model
    smp <- round(seq(dim(cmrfit$Chain)[1]/2,
                     dim(cmrfit$Chain)[1],length.out=nsamp))
    
    tmp1spec <- array(NA,c(length(cmrfit$Model$inx$inx1$specInx),6))
    tmp1nam <- c()
    for (ii in 1:length(cmrfit$Model$inx$inx1$specInx)){
      tmp1 <- cmrfit$Chain[smp,cmrfit$Model$inx$inx1$specInx[ii]];
      tmp1nam[[ii]] <- noquote(names(cmrfit$Model$inx$inx1$specInx)[[ii]]);
      if (ii==1){
        tmp1 = exp(tmp1);
        tmp1nam[[ii]] = 'Overall rate';
      }
      tmp1spec[ii,] = c(mean(tmp1),           quantile(tmp1,c(0.025,0.5,0.975)),sum(tmp1*(sign(mean(tmp1)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$inx1$specInx[ii]]))
    }
    colnames(tmp1spec) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmp1spec) <- tmp1nam;
    ### Extinction
    tmp1ext <- array(NA,c(length(cmrfit$Model$inx$inx1$extInx),6))
    tmp1nam <- c()
    for (ii in 1:length(cmrfit$Model$inx$inx1$extInx)){
      tmp1 <- cmrfit$Chain[smp,cmrfit$Model$inx$inx1$extInx[ii]];
      tmp1nam[[ii]] <- noquote(names(cmrfit$Model$inx$inx1$extInx)[[ii]]);
      if (ii==1){
        tmp1 = exp(tmp1);
        tmp1nam[[ii]] = 'Overall rate';
      }
      tmp1ext[ii,] = c(mean(tmp1),           quantile(tmp1,c(0.025,0.5,0.975)),sum(tmp1*(sign(mean(tmp1)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$inx1$extInx[ii]]))
    }
    colnames(tmp1ext) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmp1ext) <- tmp1nam;
    ## == SAMPLING
    tmp1samp <- array(NA,c(length(cmrfit$Model$inx$inx1$sampInx),6))
    tmp1nam <- c()
    for (ii in 1:length(cmrfit$Model$inx$inx1$sampInx)){
      tmp1 <- cmrfit$Chain[smp,cmrfit$Model$inx$inx1$sampInx[ii]];
      tmp1nam[[ii]] <- noquote(names(cmrfit$Model$inx$inx1$sampInx)[[ii]]);
      if (ii==1){
        tmp1 = exp(tmp1);
        tmp1nam[[ii]] = 'Overall rate';
      }
      tmp1samp[ii,] = c(mean(tmp1),           quantile(tmp1,c(0.025,0.5,0.975)),sum(tmp1*(sign(mean(tmp1)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$inx1$sampInx[ii]]))
    }
    colnames(tmp1samp) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmp1samp) <- tmp1nam;
    
    
    
    
    
    
    tmp2spec <- array(NA,c(length(cmrfit$Model$inx$inx2$specInx),6))
    tmp2nam <- c()
    for (ii in 1:length(cmrfit$Model$inx$inx2$specInx)){
      tmp2 <- cmrfit$Chain[smp,cmrfit$Model$inx$inx2$specInx[ii]];
      tmp2nam[[ii]] <- noquote(names(cmrfit$Model$inx$inx2$specInx)[[ii]]);
      if (ii==1){
        tmp2 = exp(tmp2);
        tmp2nam[[ii]] = 'Overall rate';
      }
      tmp2spec[ii,] = c(mean(tmp2),           quantile(tmp2,c(0.025,0.5,0.975)),sum(tmp2*(sign(mean(tmp2)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$inx2$specInx[ii]]))
    }
    colnames(tmp2spec) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmp2spec) <- tmp2nam;
    ### Extinction
    tmp2ext <- array(NA,c(length(cmrfit$Model$inx$inx2$extInx),6))
    tmp2nam <- c()
    for (ii in 1:length(cmrfit$Model$inx$inx2$extInx)){
      tmp2 <- cmrfit$Chain[smp,cmrfit$Model$inx$inx2$extInx[ii]];
      tmp2nam[[ii]] <- noquote(names(cmrfit$Model$inx$inx2$extInx)[[ii]]);
      if (ii==1){
        tmp2 = exp(tmp2);
        tmp2nam[[ii]] = 'Overall rate';
      }
      tmp2ext[ii,] = c(mean(tmp2),           quantile(tmp2,c(0.025,0.5,0.975)),sum(tmp2*(sign(mean(tmp2)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$inx2$extInx[ii]]))
    }
    colnames(tmp2ext) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmp2ext) <- tmp2nam;
    ## == SAMPLING
    tmp2samp <- array(NA,c(length(cmrfit$Model$inx$inx2$sampInx),6))
    tmp2nam <- c()
    for (ii in 1:length(cmrfit$Model$inx$inx2$sampInx)){
      tmp2 <- cmrfit$Chain[smp,cmrfit$Model$inx$inx2$sampInx[ii]];
      tmp2nam[[ii]] <- noquote(names(cmrfit$Model$inx$inx2$sampInx)[[ii]]);
      if (ii==1){
        tmp2 = exp(tmp2);
        tmp2nam[[ii]] = 'Overall rate';
      }
      tmp2samp[ii,] = c(mean(tmp2),           quantile(tmp2,c(0.025,0.5,0.975)),sum(tmp2*(sign(mean(tmp2)))<0)/length(smp),coda::effectiveSize(cmrfit$Chain[-c(1:(dim(cmrfit$Chain)[1]/2)),cmrfit$Model$inx$inx2$sampInx[ii]]))
    }
    colnames(tmp2samp) <- c("mean",                     "2.5%",                     "median",                     "97.5%",                     "p >/<0 ",                     "Eff SS")
    rownames(tmp2samp) <- tmp2nam;
    
    
    
    
    
    
    
    
    
    
    
    writeLines("\t \t \t=== Speciation rate parameters ===")
    writeLines("\t \t \t Clade 1")
    print(tmpspec1,digits=3)
    writeLines("\t \t \t Clade 2")
    print(tmpspec2,digits=3)
    writeLines("\t \t \t=== Extinction rate parameters ===")
    writeLines("\t \t \t Clade 1")
    print(tmpext1,digits=3)
    writeLines("\t \t \t Clade 2")
    print(tmpext2,digits=3)
    writeLines("\t \t \t=== Sampling rate parameters ===")
    writeLines("\t \t \t Clade 1")
    print(tmpsamp1,digits=3)
    writeLines("\t \t \t Clade 2")
    print(tmpsamp2,digits=3)
    
    out = list(SpecStat_Clade1= tmpspec1,
               ExtStat_Clade1 = tmpext1,
               SampStat_Clade1= tmpsamp1,
               SpecStat_Clade2= tmpspec2,
               ExtStat_Clade2 = tmpext2,
               SampStat_Clade2= tmpsamp2);
    # rownames(tmpspec) <- c('Overall rate',noquote(names(cmrfit$Model$inx$inx1$specInx)[[2]]))
    # format(tmpspec,digits=3)
    return(invisible(out));
  }
  
}
