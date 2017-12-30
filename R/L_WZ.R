#' Calculates Weaver and Zhou likelihood
#'
#' This function calculates the likelihood from Weaver and Zhou (2012)
#'
#' @details
#' The dependent variable is assumed to be normally distributed
#'
#' @param ys vector of sample values of the dependent variable
#' @param yr vector of non-sample values of the dependent variable
#' @param xs matrix of sample covariate values of dimension n by p
#' where n is the sample size (length(ys))
#' @param beta vector of regression coefficients
#' @param sd error standard deviation
#' @param cutoffs A vector of finite cutoffs defining the strata.
#'  If there are H cutoffs there are H+1 strata.
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @return likelihood or log-likelihood conditional on stratum memberships
#' @examples
#' data(population_example)
#' L.WZ(ys=sample.example$y,yr=population.example$y[population.example$s.indicators==0],
#' xs=cbind(1,sample.example$x),beta=c(1,1),sd=1,cutoffs=c(1,2,3),log=TRUE)
#' @export
L.WZ <- function(ys,yr,xs,beta,sd,cutoffs,log=F){
  H <- length(cutoffs)+1
  n <- length(ys)
  N <- n + length(yr)
  if(!is.matrix(xs)) xs <- matrix(xs,ncol=1)
  p <- ncol(xs)
  which.intervals.s <- factor(locate.points.in.intervals(ys,cutoffs)$interval.number,levels=c(1:H))
  which.intervals.r <- factor(locate.points.in.intervals(yr,cutoffs)$interval.number,levels=c(1:H))
  wt.h <- ( as.vector(table(which.intervals.r))+as.vector(table(which.intervals.s)) ) /
    as.vector(table(which.intervals.s)) / N
  wt.s <- wt.h[which.intervals.s]
  fv <- as.vector(xs%*%beta)
  l.misspecified <- sum(dnorm(x=ys,mean=fv,sd=sd,log=T))
  L.r <- rep(0,N-n)
  for(i in c(1:n))    L.r <- L.r + wt.s[i] * dnorm(x=yr,mean=fv[i],sd=sd)
  l.r.total <- sum(log(L.r))
  l <- l.misspecified + l.r.total
  if(log) return(l) else return(exp(l))
}
