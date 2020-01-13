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
#' @param cutoffs A vector of finite cutoffs defining the strata.
#'  If there are H cutoffs there are H+1 strata.
#' @param fh Sampling fraction in each stratum for the stratified part of the sample
#' @param fSRS Sampling fraction of the SRSWOR component of the sample
#' @param in.s0 A vector equal to 1 for units selected in the SRSWOR and to 0 for those
#'  in the stratified SRSWOR. Default is that there are no SRSWOR units.
#' @param type If "SL", then the sample likelihood. If "CL" then the conditional
#'             likelihood is returned, as defined by the first bracketed
#'             expression in equation (2) of Zhou et al. (2002).
#' @param estimate.var If TRUE, variance is estimated using the observed information. Defaults to FALSE
#'
#' @return a vector containing p estimated coefficients (where p=ncol(xs)) and the
#' estimated error standard deviation. If estimate.var=TRUE then variance estimates are provided.
#' @examples
#' data(population_example)
#' lm.samplik.2x(ys=sample.example$y,
#'               xs=cbind(1,sample.example$x1,sample.example$x2),
#'               cutoffs=c(0,2),fh=c(0.04,0.01,0.05))
#' @export
lm.samplik.2x <- function(ys,xs,cutoffs,fh,fSRS=0,in.s0,type="SL",estimate.var=FALSE){
  n <- length(ys)
  if(missing(in.s0)) in.s0 <- rep(0,n)
  # get starting values
  which.intervals <- locate.points.in.intervals(ys,cutoffs)
  pi.s <- fh[which.intervals$interval.number] + fSRS
  lm.piwt <- lm(ys~xs-1,weights=1/pi.s)
  sd.start <- sqrt( sum(lm.piwt$residuals^2/pi.s)/sum(1/pi.s) )
  if(!is.matrix(xs)) xs <- matrix(xs,ncol=1)
  p <- ncol(xs)
  modfit <- optim(par=c(lm.piwt$coef,log(sd.start)),
                  fn=function(par){
                    Ls.2x(ys=ys,mean=as.vector(xs%*%par[1:p]),sd=exp(par[p+1]),
                          cutoffs=cutoffs,fh=fh,fSRS=fSRS,in.s0=in.s0,type=type,
                          log=TRUE)
                  },
                  control=list(fnscale=-1) , hessian=estimate.var )
  if(!estimate.var) out <- c(modfit$par[1:p],exp(modfit$par[p+1]))
  if(estimate.var) out <- list( ests=c(modfit$par[1:p],exp(modfit$par[p+1])) ,
                                var.est=solve(-modfit$hessian) )
  out
}
