#' Calculates maximum likelihood estimators for linear regression of complex survey data
#'
#' This function fits a linear model to complex sample data by maximum likelihood
#'
#' @details
#' Add some details later.
#'
#' @param ys vector of sample values of the dependent variable
#' @param xs matrix of covariate values (number of rows must equal n where n is length of ys).
#' Intercept is not assumed, so xs should contain a column of 1s if you want an intercept in the model.
#' @param p.s A function defining the ODS sample design.
#' It must take arguments ys, yr, log and stuff, and return
#' the probability (or log of the probability) that a particular
#' sample was selected.
#' @param specs An object containing detailed specifications of the design.
#' @param pi.s n-vector (where n=length(y)) of probabilities of selection for the sample units
#' @param N the population size
#' @param R the number of simulations used in the monte carlo approximation
#'        of the likelihood
#' @return a vector containing p estimated coefficients (where p=ncol(xs)) and the
#' estimated error standard deviation. No variance estimates provided.
#' @examples
#' data(population_example)
#' lm.MLE(ys=sample.example$y,xs=cbind(1,sample.example$x),p.s=p.s.ushaped,
#'        specs=c(n0=10,tuner=0.1,return.what=1),pi.s=sample.example$pi,N=100,R=10)
#' @export
lm.MLE <- function(ys,xs,p.s,specs,pi.s,N,R){
  # get starting values
  lm.piwt <- lm(ys~xs-1,weights=1/pi.s)
  sd.start <- sqrt( sum(lm.piwt$residuals^2/pi.s)/sum(1/pi.s) )
  if(!is.matrix(xs)) xs <- matrix(xs,ncol=1)
  p <- ncol(xs)
  n <- length(ys)
  all.resamp.x <- matrix(integer(length=(N-n)*R),nrow=R,ncol=N-n)
  for(r in c(1:R)) all.resamp.x[r,] <- sample(n,N-n,replace=TRUE)
  all.errors.y <- matrix(data=rnorm(n=(N-n)*R),nrow=R,ncol=N-n)
  modfit <- optim(par=c(lm.piwt$coef,log(sd.start)),
                  fn=function(par,ys,xs,N,p.s,specs,pi.s,R,all.resamp.x,all.errors.y){
                    L(beta=par[1:p],sd=exp(par[p+1]),ys=ys,xs=xs,N=N,p.s=p.s,
                      specs=specs,log=TRUE,pi.s=pi.s,R=R,all.resamp.x=all.resamp.x,all.errors.y = all.errors.y)
                  },
                  ys=ys,xs=xs,N=N,p.s=p.s,specs=specs,pi.s=pi.s,R=R,
                  control=list(fnscale=-1) )
  c(modfit$par[1:p],exp(modfit$par[p+1]))
}
