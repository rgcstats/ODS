#' Calculates maximum sample likelihood estimators for linear regression
#'
#' This function fits a linear model to complex sample data by maximum sample likelihood
#'
#' @details
#' See Krieger and Pfefferman's paper.
#' The dependent variable is assumed to be normally distributed.
#'
#' @param ys vector of sample values of the dependent variable
#' @param xs matrix of covariate values (number of rows must equal n where n is length of ys).
#' Intercept is not assumed, so xs should contain a column of 1s if you want an intercept in the model.
#' @param pi.fn a function taking a single argument y and returning the probability
#' of selection
#' @param pi.h selection probability in each stratum (H-vector)
#' @param cutoffs vector of H-1 cutoffs on Y defining the strata
#' @param estimate.var If TRUE, variance is estimated using the observed information. Defaults to FALSE
#'
#' @return a vector containing p estimated coefficients (where p=ncol(xs)) and the
#' estimated error standard deviation. If estimate.var=TRUE then variance estimates are provided.
#' @examples
#' data(population_example)
#' lm.samplik(ys=sample.example$y,xs=cbind(1,sample.example$x),
#' pi.fn=function(y,log=FALSE){
#' out <- sqrt(0.1+0.9*(y-5)^2/1)
#' if(log) return(log(out))
#' if(!log) return(out)})
#' @export
lm.samplik.2x.poisson.strat <- function(ys,xs,pi.fn,pi.h,cutoffs,estimate.var=FALSE){
  # get starting values
  pi.s <- pi.fn(ys)
  lm.piwt <- lm(ys~xs-1,weights=1/pi.s)
  sd.start <- sqrt( sum(lm.piwt$residuals^2/pi.s)/sum(1/pi.s) )
  if(!is.matrix(xs)) xs <- matrix(xs,ncol=1)
  p <- ncol(xs)
  modfit <- optim(par=c(lm.piwt$coef,log(sd.start)),
                  fn=function(par){
                    Ls.2x.poisson.strat(ys=ys,mean=as.vector(xs%*%par[1:p]),sd=exp(par[p+1]),
                       pi.h=pi.h,cutoffs=cutoffs,log=TRUE)
                  },
                  control=list(fnscale=-1),hessian=estimate.var )
  if(!estimate.var) out <- c(modfit$par[1:p],exp(modfit$par[p+1]))
  if(estimate.var) out <- list( ests=c(modfit$par[1:p],exp(modfit$par[p+1])) ,
                                var.est=solve(-modfit$hessian) )
  out
}
