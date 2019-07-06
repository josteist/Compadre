#' @export
# make a 'summary' function for a fit.

summary.CMR_fit <- function(cmrfit,nsamp = 1e4){
  # quick stats function for a named variable.
  # Output mean,p2.5,p50,p97.5,p>/<(cmrfit$Chain)[1]/2,

  smp <- round(seq(dim(cmrfit$Chain)[1]/2,
             dim(cmrfit$Chain)[1],length.out=nsamp))
  writeLines("\t \t \t=== Speciation rate parameters ===")
  writeLines(c("       ","\t\t\t",
               "mean",
               "2.5%",
               "median",
               "97.5%",
               "p >/<0 ",
               "Eff SS"),sep="\t")
  cat("\n")
  for (ii in 1:length(cmrfit$Model$inx$specInx)){
    tmp <- cmrfit$Chain[smp,cmrfit$Model$inx$specInx[ii]];
    tmpnam <- noquote(names(cmrfit$Model$inx$specInx)[[ii]]);
    if (ii==1){
      tmp = exp(tmp);
      tmpnam = 'Overall rate';

    }

    writeLines(c(tmpnam,
                 ifelse(nchar(names(cmrfit$Model$inx$specInx)[[ii]])>7,"","\t"),
                 ifelse(nchar(names(cmrfit$Model$inx$specInx)[[ii]])>14,"","\t"),
                 format(c(mean(tmp),
                          quantile(tmp,c(0.025,0.5,0.975)),
                          sum(tmp*(sign(mean(tmp)))<0)/length(smp),
                          coda::effectiveSize(cmrfit$Chain[smp,cmrfit$Model$inx$specInx[ii]])),digits=2,nsmall=2)),sep="\t")
    cat("\n")
  }
  cat("\n")
  writeLines("\t \t \t=== Extinction rate parameters ===")
  writeLines(c("       ","\t\t\t",
               "mean",
               "2.5%",
               "median",
               "97.5%",
               "p >/<0 ",
               "Eff SS"),sep="\t")
  cat("\n")
  # Should we 'unlog' intercept to get the mean 'rate'?
  for (ii in 1:length(cmrfit$Model$inx$extInx)){
    tmp <- cmrfit$Chain[smp,cmrfit$Model$inx$extInx[ii]];
    tmpnam <- noquote(names(cmrfit$Model$inx$extInx)[[ii]]);
    if (ii==1){
      tmp = exp(tmp);
      tmpnam = 'Overall rate';

    }
    writeLines(c(tmpnam,
                 ifelse(nchar(names(cmrfit$Model$inx$extInx)[[ii]])>7,"","\t"),
                 ifelse(nchar(names(cmrfit$Model$inx$extInx)[[ii]])>14,"","\t"),
                 format(c(mean(tmp),
                          quantile(tmp,c(0.025,0.5,0.975)),
                          sum(tmp*(sign(mean(tmp)))<0)/length(smp),
                          coda::effectiveSize(cmrfit$Chain[smp,cmrfit$Model$inx$extInx[ii]])),digits=2,nsmall=2)),sep="\t")
    cat("\n")
  }
  cat("\n")
  writeLines("\t \t \t=== Sampling rate parameters ===")
  writeLines(c("       ","\t\t\t",
               "mean",
               "2.5%",
               "median",
               "97.5%",
               "p >/<0 ",
               "Eff SS"),sep="\t")
  cat("\n")
  for (ii in 1:length(cmrfit$Model$inx$sampInx)){
    tmp <- cmrfit$Chain[smp,cmrfit$Model$inx$sampInx[ii]];
    tmpnam <- noquote(names(cmrfit$Model$inx$sampInx)[[ii]]);
    if (ii==1){
      tmp = exp(tmp);
      tmpnam = 'Overall rate';

    }

    writeLines(c(tmpnam,
                 ifelse(nchar(names(cmrfit$Model$inx$sampInx)[[ii]])>7,"","\t"),
                 ifelse(nchar(names(cmrfit$Model$inx$sampInx)[[ii]])>14,"","\t"),
                 format(c(mean(tmp),
                          quantile(tmp,c(0.025,0.5,0.975)),
                          sum(tmp*(sign(mean(tmp)))<0)/length(smp),
                          coda::effectiveSize(cmrfit$Chain[smp,cmrfit$Model$inx$sampInx[ii]])),digits=2,nsmall=2)),sep="\t")
    cat("\n")
  }
}

# mysf(fa)
# THis could work (the length of stats output is increased by 1 if neg)
# have a switch for length of 'driver name' if > 6 chars
