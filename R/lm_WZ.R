#' Linear model fitting for stratified survey data by maximising the estimated likelihood of yu and xs
#'
#' This function fits a linear model to survey data stratified by
#' the dependent variable, by maximising the estimated likelihood
#' from Weaver and Zhou (2012).
#'
#' @param ys vector of sample values of the dependent variable
#' @param yr vector of non-sample values of the dependent variable
#' @param xs matrix of covariate values (number of rows must equal n where n is length of ys).
#' Intercept is not assumed, so xs should contain a column of 1s if you want an intercept in the model.
#' @param cutoffs A vector of finite cutoffs defining the strata.
#'  If there are H cutoffs there are H+1 strata.
#' @param pi.s n-vector (where n=length(y)) of probabilities of selection for the sample units
#' @return a vector containing p estimated coefficients (where p=ncol(xs)) and the
#' estimated error standard deviation. No variance estimates provided.
#' @examples
#' data(population_example)
#' lm.WZ(ys=sample.example$y,yr=population.example$y[population.example$s.indicators==0],
#' xs=cbind(1,sample.example$x),cutoffs=c(1,2,3),pi.s=sample.example$pi)
#' @export
lm.WZ <- function(ys,yr,xs,cutoffs,pi.s){
  n <- length(ys)
  if(missing(pi.s)) pi.s <- rep(1,n)
  lm.piwt <- lm(ys~xs-1,weights=1/pi.s)
  sd.start <- sqrt( sum(lm.piwt$residuals^2/pi.s)/sum(1/pi.s) )
  if(!is.matrix(xs)) xs <- matrix(xs,ncol=1)
  p <- ncol(xs)
  modfit <- optim(par=c(lm.piwt$coef,log(sd.start)),
                  fn=function(par,ys,yr,xs,cutoffs){
                    L.WZ(ys=ys,yr=yr,xs=xs,beta=par[1:p],sd=exp(par[p+1]),
                         cutoffs=cutoffs,log=T)
                  },
                  ys=ys,yr=yr,xs=xs,cutoffs=cutoffs,control=list(fnscale=-1) )
  c(modfit$par[1:p],exp(modfit$par[p+1]))
}
