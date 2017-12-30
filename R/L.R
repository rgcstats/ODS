#' Calculates the estimated likelihood for a regression model under a general sample design
#'
#' This function calculates the likelihood
#'
#' @details
#' Add some details later.
#'
#' @param beta vector of regression coefficients
#' @param sd error standard deviation
#' @param ys vector of sample values of the dependent variable
#' @param xs matrix of sample covariate values of dimension n by p
#' where n is the sample size (length(ys))
#' @param N the population size
#' @param p.s A function defining the ODS sample design.
#' It must take arguments ys, yr, log and stuff, and return
#' the probability (or log of the probability) that a particular
#' sample was selected.
#' @param specs An object containing detailed specifications of the design.
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @param pi.s the probabilities of selection of the sample units
#' @param R the number of simulations used in the monte carlo approximation
#'        of the likelihood
#' @param all.resamp.x Optional. An R by N-n integer-valued matrix
#'        whose rows are SRSWR resamples from 1:n
#' @param all.errors.y Optional. An R by N-n real-valued matrix
#' of standard normal realisations.
#' @param verbose If TRUE, the function outputs information to the console.
#' @return The likelihood or log-likelihood.
#' @examples
#' data(population_example)
#' L(beta=c(4,1),N=100,sd=0.8,ys=sample.example$y,xs=cbind(1,sample.example$y),
#' p.s=p.s.ushaped,specs=c(n0=10,tuner=0.1,return.what=1),log=TRUE,
#' pi.s=sample.example$pi)
#' @export
L <- function(beta,sd,ys,xs,N,p.s,specs=NULL,log=FALSE,pi.s,R=100,
              all.resamp.x,all.errors.y,verbose=FALSE){
  if(!is.matrix(xs)) xs <- matrix(xs,ncol=1)
  p <- ncol(xs)
  n <- length(ys)
  l.misspecified <- sum( dnorm(x=ys,mean=as.vector(xs%*%beta),sd=sd,log=T) )
  all.log.p.s <- rep(NA,R)
  for(r in c(1:R)){
    if(missing(all.resamp.x)) resamp.x <- sample(x=1:n,size=N-n,prob=1/pi.s,replace=T)
    if(!missing(all.resamp.x)) resamp.x <- all.resamp.x[r,]
    if(missing(all.errors.y)) errors.y <- rnorm(n=N-n)
    if(!missing(all.errors.y)) errors.y <- all.errors.y[r,]
    xr.sim <- xs[resamp.x,]
    yr.sim <- as.vector(xr.sim %*% beta) + errors.y*sd
    all.log.p.s[r] <- p.s(ys=ys,yr=yr.sim,specs=specs,log=T)
  }
  K <- median(all.log.p.s)
  log.sampdesign.factor <- log( mean(exp(all.log.p.s-K)) ) + K
  l <- l.misspecified + log.sampdesign.factor
  if(verbose) cat("l.misspecified=",l.misspecified," ; l=",l," beta=(",paste(beta,sep=""),") ; sd=",sd,"\n")
  if(log) return(l)
  if(!log) return(exp(l))
}
