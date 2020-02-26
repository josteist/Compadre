#' A function to generate sensible initial values for MCMC estimation.
#'
#' getInits generates initial values for a complex CMR model. Essentially the global parameters (mean rates) are estimated optimizing the probability density function, assuming all random effects are 0. Then each random effect is estimated independently, assuming all other random effects are 0. This is not proper optimization, but work reasonably well for initial values.
#'
#' @param cmrModel A CMR model structure with 1 clade.
#' @return The function returns relatively decent initial values for MCMC estimation.
#' @export
#' @examples mod1 <- make.BayesCMR(Obs)
#' x0 <- getInit(mod1)
#' fit <- MCMC_CMR(mod1,x0=x0)




getInit <- function(cmrModel){
  # THis function will generate probable initial values for the MCMC chain.
  # The main parameters will be estimated by optim ignoring the random effects afterwhich all the RE will be optimzed independently (1 dimension at a time).
  # Change to also include drivers, but NOT variances.
  if (length(cmrModel$inx)==2){
    # it's a 2-clade omdel
    x_1 <- c(optim(c(-1,-1.1,-1.2),function(x){-cmrModel$Clade1Mod$probfun(c(x[1],x[2],x[3],rep(0,cmrModel$Clade1Mod$npar-3)))})$par,
             rep(0,cmrModel$Clade1Mod$npar-3))
    xix <- sample(unlist(cmrModel$Clade1Mod$inx[5:7]))
    if (length(xix)>0){
      xres <- sapply(xix,function(jj){
        optim(0,function(yy){-cmrModel$Clade1Mod$probfun(c(x_1[1:(jj-1)],yy,x_1[(jj+1):cmrModel$Clade1Mod$npar]))},
              method='Brent',lower=-10,upper=10)$par})
      x_1[xix]=xres;
    }

    x_2 <- c(optim(c(-1,-1.1,-1.2),function(x){-cmrModel$Clade2Mod$probfun(c(x[1],x[2],x[3],rep(0,cmrModel$Clade2Mod$npar-3)))})$par,
             rep(0,cmrModel$Clade2Mod$npar-3))
    xix <- sample(unlist(cmrModel$Clade2Mod$inx[5:7]))
    if (length(xix)>0){
      xres <- sapply(xix,function(jj){
        optim(0,function(yy){-cmrModel$Clade2Mod$probfun(c(x_2[1:(jj-1)],yy,x_2[(jj+1):cmrModel$Clade2Mod$npar]))},
              method='Brent',lower=-10,upper=10)$par})
      x_2[xix]=xres;
    }
    x = c(x_1,x_2);


  } else {

    if (cmrModel$npar>3){
      x <- c(optim(c(-1,-1.1,-1.2),function(x){-cmrModel$probfun(c(x[1],x[2],x[3],rep(0,cmrModel$npar-3)))})$par,
             rep(0,cmrModel$npar-3))
      xix <- sample(unlist(cmrModel$inx[5:7]))
      if (length(xix)>0){
        xres <- sapply(xix,function(jj){
          optim(0,function(yy){-cmrModel$probfun(c(x[1:(jj-1)],yy,x[(jj+1):cmrModel$npar]))},
                method='Brent',lower=-10,upper=10)$par})
        x[xix]=xres;
      }
    } else {
      x <- optim(c(-1,-1.1,-1.2),function(x){-cmrModel$probfun(x)})$par
    }
  }
  return(x0=x)
}
