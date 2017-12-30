#' Calculates maximum likelihood estimators for linear regression of complex survey data
#'
#' This function fits a linear model to complex sample data by maximum likelihood
#'
#' @details
#' Add some details later.
#'
#' @param ys vector of sample values of the dependent variable
#' @param x1s sample values for covariate 1 (continuous)
#' @param x2s sample values for covariate 2 (binary)
#' @param pi.fn The probability that a unit is selected as a function of y.
#' @param pi.s n-vector (where n=length(y)) of probabilities of selection for the sample units
#' @param Rquad number of quadrature points used in calculating the likelihood
#' @return a vector containing p estimated coefficients (where p=ncol(xs)) and the
#' estimated error standard deviation. No variance estimates provided.
#' @examples
#' data(population_example)
#' eg.ybar.pi <- sum(sample.example$y/sample.example$pi) / sum(1/sample.example$pi)
#' eg.s2.pi <- sum((sample.example$y-eg.ybar.pi)^2/sample.example$pi)/sum(1/sample.example$pi)
#' eg.pi.fn <- function(y){ out <- sqrt( 0.1 + 0.9*(y-eg.ybar.pi)^2/eg.s2.pi )
#' out <- out * mean(sample.example$pi) /
#'        mean(sqrt( 0.1 + 0.9*(sample.example$y-eg.ybar.pi)^2/eg.s2.pi ))
#' out  }
#' lm.MLE.2x.poisson(ys=sample.example$y,x1s=sample.example$x1,x2s=sample.example$x2,pi.fn=eg.pi.fn,
#'        pi.s=sample.example$pi,R=10)
#' @export
lm.MLE.2x.poisson.strat <- function(ys,x1s,x2s,pi.s,pi.h,cutoffs){
  # get starting values
  lm.piwt <- lm(ys~x1s+x2s,weights=1/pi.s)
  sd.start <- sqrt( sum(lm.piwt$residuals^2/pi.s)/sum(1/pi.s) )
  mu0.x1.start <- sum(x1s*(1-x2s)/pi.s) / sum((1-x2s)/pi.s)
  mu1.x1.start <- sum(x1s*x2s/pi.s) / sum(x2s/pi.s)
  sd.x1.start <- sqrt( sum( (x1s-(1-x2s)*mu0.x1.start-x2s*mu1.x1.start)^2 / pi.s ) / sum(1/pi.s) )
  prob.x2.start <- sum(x2s/pi.s)/sum(1/pi.s)
  n <- length(ys)
  # optimise over parameters: beta (3 vector), log(sd), logit(prob.x2) , mu0.x1, mu1.x1, log(sd.x1),
  startpar <- c( lm.piwt$coef , log(sd.start), log(prob.x2.start/(1-prob.x2.start)) ,
                  mu0.x1.start , mu1.x1.start , log(sd.x1.start) )
  modfit <- optim(par=startpar,
                  fn=function(par){
                    L.2x.poisson.strat(beta=par[1:3],sd=exp(par[4]),prob.x2=exp(par[5])/(1+exp(par[5])),
                         mu0.x1=par[6] , mu1.x1=par[7] ,sd.x1=exp(par[8]),
                         ys=ys,x1s=x1s,x2s=x2s,pi.fn=pi.fn,
                         log=TRUE,pi.h=pi.h,cutoffs=cutoffs)
                  },
                  control=list(fnscale=-1,maxit=1000) )
  c(modfit$par[1:3],exp(modfit$par[4]),
    exp(modfit$par[5])/(1+exp(modfit$par[5])) , modfit$par[6:7] ,
    exp(modfit$par[8]) )
}
