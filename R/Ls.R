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
#' @return a quantity proportional to the density of ys (or the log thereof)
#' @examples
#' data(population_example)
#' Ls(ys=sample.example$y,mean=5+0.5*sample.example$x,sd=1,log=TRUE,
#' pi.fn=function(y,log=FALSE){
#' out <- sqrt(0.1+0.9*(y-5)^2/1)
#' if(log) return(log(out))
#' if(!log) return(out)})
#' @export
Ls <- function(ys,mean,sd,pi.fn,log=FALSE,R=100){
  n <- length(ys)
  if(length(sd)==1) sd <- rep(sd,n)
  zsim <- qnorm(p=c(1:R)/(R+1))
  ylong <- rep(mean,each=R) + rep(sd,each=R)*rep(zsim,times=n)
  integrand <- pi.fn(ylong)
  integral <- tapply(integrand,rep(c(1:n),each=R),mean)
  ls <- sum( dnorm(ys,mean,sd,log=T) + log(pi.fn(ys)) - log(integral) )
  if(log) return(ls) else return(exp(ls))
}

