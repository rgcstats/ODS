#' Calculates the profiled empirical likelihood
#'
#' This function calculates the profiled empirical likelihood (or rather something proportional to it)
#'
#' @details
#' The dependent variable is assumed to be normally distributed
#'
#' @param ys vector of sample values of the dependent variable
#' @param mean vector of means of dependent variable
#' @param sd vector of standard deviations of dependent variable
#' @param pi.fn a function taking a single argument y and returning the probability
#' of selection
#' @param N the population size
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @param R number of points used in quadrature calculation of integral
#' @param tol a factor controlling error-tolerance in integration (defaults to 1)
#'
#' @return a quantity proportional to the density of ys (or the log thereof)
#'
#' @import stats
#' @import fastGHQuad
#'
#' @examples
#' data(population_example)
#' L.el(ys=sample.example$y,mean=5+0.5*sample.example$x,sd=1,log=TRUE,N=1000,
#' pi.fn=function(y,log=FALSE){
#' out <- pmin(1,0.05*sqrt(0.1+0.9*(y-5)^2/1))
#' if(log) return(log(out))
#' if(!log) return(out)})
#' @export
L.el <- function(ys,mean,sd,pi.fn,N,log=FALSE,R=100,tol=1){
  n <- length(ys)
  if(length(sd)==1) sd <- rep(sd,n)
  # first get empirical likelihood estimators of p.i.x (mean is assumed to be E[Y|x])
  #zsim <- qnorm(p=c(1:R)/(R+1))
  #ylong <- rep(mean,each=R) + rep(sd,each=R)*rep(zsim,times=n)
  #integrand <- pi.fn(ylong)
  #mean.pi.given.x <- tapply(integrand,rep(c(1:n),each=R),mean)
  rule <- fastGHQuad::gaussHermiteData(R)
  zsim <- rule$x
  ylong <- rep(mean,each=R) + rep(sd,each=R)*sqrt(2)*rep(zsim,times=n)
  wlong <- rep(rule$w,times=n)/sqrt(pi)
  integrand <- pi.fn(ylong)*wlong
  mean.pi.given.x <- tapply(integrand,rep(c(1:n),each=R),sum)
  #mean.pi.given.x <- rep(NA,n)
  #for(i in c(1:n)){
  #  mean.pi.given.x[i] <- integrate(f=function(y){dnorm(y,mean=mean[i],sd=sd[i])*pi.fn(y)},
  #                                  lower=-Inf,upper=Inf,rel.tol=tol*.Machine$double.eps^0.25)$value
  #}
  m.i <- 1 - mean.pi.given.x
  if(max(m.i)==min(m.i)) lambda<-max(m.i) else{
    objfn <- function(lambda,m.i,N,n){
      p.i <- 1 / ( N - (N-n)*m.i/lambda )
      ((lambda-sum(m.i*p.i))/lambda)^2
    }
    lowerlim <- pmax( min(m.i) , N / (N-1) * (1-n/N) * max(m.i) )
    lowerlim <- pmin(lowerlim+0.005,max(m.i)-0.001)
    lambda <- optimise(f=objfn,interval=c(lowerlim,max(m.i)),m.i=m.i,N=N,n=n)$minimum
  }
  p.i.x <- 1 / ( N - (N-n)*m.i/lambda )
  p.i.x <- pmax(1e-6,p.i.x)
  p.i.x <- p.i.x / sum(p.i.x)
  l <- sum( dnorm(ys,mean,sd,log=T) ) + sum(log(p.i.x)) + (N-n)*log(sum(p.i.x*m.i))
  if(log) return(l) else return(exp(l))
}

