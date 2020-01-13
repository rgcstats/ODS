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
#' @param p.s A function defining the ODS sample design.
#' It must take arguments ys, yr, log and stuff, and return
#' the probability (or log of the probability) that a particular
#' sample was selected.
#' @param pi.fn Optional. The approximate probability that a unit is selected as a function of y. If provided,
#'              a more efficient algorithm is used.
#' @param specs An object containing detailed specifications of the design.
#' @param pi.s n-vector (where n=length(y)) of probabilities of selection for the sample units
#' @param N the population size
#' @param R the number of simulations used in the monte carlo approximation
#'        of the likelihood
#' @param Rquad number of quadrature points
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
#' lm.MLE.2x(ys=sample.example$y,x1s=sample.example$x1,x2s=sample.example$x2,p.s=p.s.ushaped,
#'        specs=c(n0=10,tuner=0.1,return.what=1),pi.s=sample.example$pi,N=1000,R=10)
#' @export
lm.EM.2x <- function(ys,x1s,x2s,p.s,specs,pi.s,pi.fn,N,R,Rquad){
  # get starting values
  lm.piwt <- lm(ys~x1s+x2s,weights=1/pi.s)
  sd.start <- sqrt( sum(lm.piwt$residuals^2/pi.s)/sum(1/pi.s) )
  mu0.x1.start <- sum(x1s*(1-x2s)/pi.s) / sum((1-x2s)/pi.s)
  mu1.x1.start <- sum(x1s*x2s/pi.s) / sum(x2s/pi.s)
  sd.x1.start <- sqrt( sum( (x1s-(1-x2s)*mu0.x1.start-x2s*mu1.x1.start)^2 / pi.s ) / sum(1/pi.s) )
  prob.x2.start <- sum(x2s/pi.s)/sum(1/pi.s)
  n <- length(ys)
  if(is.null(pi.fn)) expand <- 1 else expand <- 4
  all.errors.y <- matrix(data=rnorm(n=(N-n)*R*expand),nrow=R,ncol=(N-n)*expand)
  all.errors.x1 <- matrix(data=rnorm(n=(N-n)*R*expand),nrow=R,ncol=(N-n)*expand)
  all.unif.x2 <- matrix(data=runif(n=(N-n)*R*expand),nrow=R,ncol=(N-n)*expand)
  all.subsamp.unif <- matrix(data=runif(n=(N-n)*R*expand),nrow=R,ncol=(N-n)*expand)
  # optimise over parameters: beta (3 vector), log(sd), logit(prob.x2) , mu0.x1, mu1.x1, log(sd.x1),
  startpar0 <- c( lm.piwt$coef , log(sd.start), log(prob.x2.start/(1-prob.x2.start)) ,
                  mu0.x1.start , mu1.x1.start , log(sd.x1.start) )
  if(is.null(startpar)) startpar <- startpar0
  if(any(is.na(startpar))) startpar[is.na(startpar)] <- startpar0[is.na(startpar)]
  l.start <- L.2x(beta=startpar[1:3],sd=exp(startpar[4]),prob.x2=exp(startpar[5])/(1+exp(startpar[5])),
                  mu0.x1=startpar[6] , mu1.x1=startpar[7] ,sd.x1=exp(startpar[8]),
                  ys=ys,x1s=x1s,x2s=x2s,N=N,p.s=p.s,pi.fn=pi.fn,
                  specs=specs,log=TRUE,pi.s=pi.s,R=R,
                  all.errors.y = all.errors.y,all.errors.x1=all.errors.x1,
                  all.unif.x2=all.unif.x2,
                  all.subsamp.unif=all.subsamp.unif)
  modfit <- optim(par=startpar,
                  fn=function(par,ys,x1s,x2s,N,p.s,pi.fn,specs,pi.s,R,all.errors.y,
                              all.errors.x1,all.unif.x2,all.subsamp.unif){
                    L.2x(beta=par[1:3],sd=exp(par[4]),prob.x2=exp(par[5])/(1+exp(par[5])),
                         mu0.x1=par[6] , mu1.x1=par[7] ,sd.x1=exp(par[8]),
                         ys=ys,x1s=x1s,x2s=x2s,N=N,p.s=p.s,pi.fn=pi.fn,
                         specs=specs,log=TRUE,pi.s=pi.s,R=R,
                         all.errors.y = all.errors.y,all.errors.x1=all.errors.x1,
                         all.unif.x2=all.unif.x2,
                         all.subsamp.unif=all.subsamp.unif)
                  },
                  ys=ys,x1s=x1s,x2s=x2s,N=N,p.s=p.s,pi.fn=pi.fn,
                  specs=specs,pi.s=pi.s,R=R,
                  all.errors.y = all.errors.y,all.errors.x1=all.errors.x1,
                  all.unif.x2=all.unif.x2,all.subsamp.unif=all.subsamp.unif,
                  control=list(fnscale=-1,maxit=1000) )
  c(modfit$par[1:3],exp(modfit$par[4]),
    exp(modfit$par[5])/(1+exp(modfit$par[5])) , modfit$par[6:7] ,
    exp(modfit$par[8]) )
}
