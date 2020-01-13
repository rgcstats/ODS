#' Calculates the sample likelihood
#'
#' This function calculates the sample likelihood (or rather something proportional to it)
#'
#' @details
#' See Krieger and Pfefferman's paper.
#' The dependent variable is assumed to be normally distributed
#'
#' @param ys vector of sample values of the dependent variable
#' @param mean vector of means of dependent variable
#' @param sd vector of standard deviations of dependent variable
#' @param pi.fn a function taking a single argument y and returning the probability
#' of selection
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @param R number of points used in quadrature calculation of integral
#' @param method Character. "GH" for Gauss-Hermite quadrature, "q" for quantile-based integration,
#' "i" to use integrate function (slow).
#' @param subdivisions Number of subdivisions to use in integration when method="i"
#'
#' @return a quantity proportional to the density of ys (or the log thereof)
#'
#' @import fastGHQuad
#'
#' @examples
#' data(population_example)
#' Ls(ys=sample.example$y,mean=5+0.5*sample.example$x,sd=1,log=TRUE,
#' pi.fn=function(y,log=FALSE){
#' out <- sqrt(0.1+0.9*(y-5)^2/1)
#' if(log) return(log(out))
#' if(!log) return(out)})
#' @export
Ls <- function(ys,mean,sd,pi.fn,log=FALSE,R=150,method="GH",subdivisions=100){
  n <- length(ys)
  if(length(sd)==1) sd <- rep(sd,n)
  if(method=="q"){
    zsim <- qnorm(p=c(1:R)/(R+1))
    ylong <- rep(mean,each=R) + rep(sd,each=R)*rep(zsim,times=n)
    integrand <- pi.fn(ylong)
    mean.pi.given.x <- tapply(integrand,rep(c(1:n),each=R),mean)
  }
  if(method=="GH"){
    rule <- fastGHQuad::gaussHermiteData(R)
    zsim <- rule$x
    ylong <- rep(mean,each=R) + rep(sd,each=R)*sqrt(2)*rep(zsim,times=n)
    wlong <- rep(rule$w,times=n)/sqrt(pi)
    integrand <- pi.fn(ylong)*wlong
    mean.pi.given.x <- tapply(integrand,rep(c(1:n),each=R),sum)
  }
  if(method=="i"){
    integrand.fn <- function(z,mean,sd){dnorm(z)*pi.fn(mean+sd*z)}
    mean.pi.given.x <- rep(NA,n)
    for(j in c(1:n)){
      mean.pi.given.x[j] <- integrate(f=integrand.fn,lower=-Inf,upper=Inf,
                                      mean=mean[j],sd=sd[j],
                                      subdivisions=subdivisions)$value
    }
  }
  ls <- sum( dnorm(ys,mean,sd,log=T) + log(pi.fn(ys)) - log(mean.pi.given.x) )
  if(log) return(ls) else return(exp(ls))
}

