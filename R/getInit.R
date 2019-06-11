#' A function to generate sensible initial values for MCMC estamation.
#'
#' getInits generates initial values for a complex CMR model. Essentially the global parameters (mean rates) are estimated optimizing the probability density function, assuming all random effects are 0. Then each random effect is estimated independently, assuming all other random effects are 0. This is not proper optimization, but work reasonably well for initial values.
#' @param cmrModel A CMR model structure.
#' @return The function returns relatively decent initial values for MCMC estimation.
#' @export
#' @examples mod1 <- make.BayesCMR(Obs)
#' x0 <- getInit(mod1)
#' fit <- MCMC_CMR(mod1,x0=x0)



# Getting initial values?

getInit <- function(cmrModel){
  # THis function will generate probable initial values for the MCMC chain.
  # The main parameters will be estimated by optim ignoring the random effects afterwhich all the RE will be optimzed independently (1 dimension at a time).
  # This will not be proper optimization, but be better initial values for analysis.
  # Change to also include drivers, but NOT variances.
  x <- c(optim(c(-1,-1.1,-1.2),function(x){-cmrModel$probfun(c(x[1],x[2],x[3],rep(0,cmrModel$npar-3)))})$par,
         1,1,1,rep(0,cmrModel$npar-6))
  xres <- sapply((1+cmrModel$nhps):(cmrModel$npar),function(jj){
    optim(0,function(yy){-cmrModel$probfun(c(x[1:(jj-1)],yy,x[(jj+1):cmrModel$npar]))},
          method='Brent',lower=-10,upper=10)$par})
  x[7:cmrModel$npar]=xres;
  return(x0=x)
}

#
# x <- fx1b$Chain[1,]
# x2 <-   sapply(7:mx1$npar,
#                function(jj){optim(c(-1),function(y){-mx1$probfun(c(x[1:jj],y,x[(jj+1):mx1$npar]))},method='Brent',lower=-10,upper=10)$par})
#
