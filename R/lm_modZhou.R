#' Linear model fitting for stratified survey data by maximum conditional likelihood
#'
#' This function fits a linear model to survey data stratified by
#' the dependent variable, by maximising a conditional likelihood
#' which is a modification of the one in Zhou et al. (2002).
#'
#' @param ys vector of sample values of the dependent variable
#' @param xs matrix of covariate values (number of rows must equal n where n is length of ys).
#' Intercept is not assumed, so xs should contain a column of 1s if you want an intercept in the model.
#' @param cutoffs A vector of finite cutoffs defining the strata.
#'  If there are H cutoffs there are H+1 strata.
#' @param pi.s n-vector (where n=length(y)) of probabilities of selection for the sample units
#' @param in.s0 A vector equal to 1 for units selected in the SRSWOR and to 0 for those
#'  in the stratified SRSWOR. Default is that there are no SRSWOR units.
#' @return a vector containing p estimated coefficients (where p=ncol(xs)) and the
#' estimated error standard deviation. No variance estimates provided.
#' @examples
#' data(population_example)
#' lm.modZhou(ys=sample.example$y,xs=cbind(1,sample.example$x),
#' cutoffs=c(0,2),pi.s=sample.example$pi)
#' @export
lm.modZhou <- function(ys,xs,cutoffs,pi.s,in.s0){
  # cutoffs must be a vector of finite stratum cutoffs
  #  if there are H strata there should be (H-1) cutoffs
  # pi.s used only to get starting values
  n <- length(ys)
  if(missing(pi.s)) pi.s <- rep(1,n)
  if(missing(in.s0)) in.s0 <- rep(0,n)
  lm.piwt <- lm(ys~xs-1,weights=1/pi.s)
  sd.start <- sqrt( sum(lm.piwt$residuals^2/pi.s)/sum(1/pi.s) )
  if(!is.matrix(xs)) xs <- matrix(xs,ncol=1)
  p <- ncol(xs)
  modfit <- optim(par=c(lm.piwt$coef,log(sd.start)),
                  fn=function(par,ys,cutoffs){
                    L.modZhou(ys=ys,mean=as.vector(xs%*%par[1:p]),sd=exp(par[p+1]),
                              cutoffs=cutoffs,log=TRUE,in.s0=in.s0)
                  },
                  ys=ys,cutoffs=cutoffs,control=list(fnscale=-1) )
  c(modfit$par[1:p],exp(modfit$par[p+1]))
}
