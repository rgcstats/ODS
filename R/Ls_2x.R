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
#' @param cutoffs A vector of finite cutoffs defining the strata.
#'  If there are H cutoffs there are H+1 strata.
#' @param fh Sampling fraction in each stratum for the stratified part of the sample
#' @param fSRS Sampling fraction of the SRSWOR component of the sample
#' @param in.s0 A vector equal to 1 for units selected in the SRSWOR and to 0 for those
#'  in the stratified SRSWOR. Default is that there are no SRSWOR units.
#' @param type If "SL", then the sample likelihood. If "CL" then the conditional
#'             likelihood is returned, as defined by the first bracketed
#'             expression in equation (2) of Zhou et al. (2002).
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#'
#' @return a quantity proportional to the density of ys (or the log thereof)
#'
#' @examples
#' data(population_example)
#' Ls.2x(mean=1+0.5*sample.example$x1-0.5*sample.example$x2,
#'       sd=1,ys=sample.example$y,cutoffs=c(0,2),fh=c(0.04,0.01,0.05),log=TRUE)
#' @export
Ls.2x <- function(ys,mean,sd,cutoffs,fh,fSRS=0,in.s0,type="SL",log=FALSE){
  n <- length(ys)
  if(length(sd)==1) sd <- rep(sd,n)
  H <- length(cutoffs)+1
  if(missing(in.s0)) in.s0 <- rep(0,n)
  in.s1 <- 1 - in.s0
  n0 <- sum(in.s0)
  n1 <- sum(in.s1)
  which.intervals <- locate.points.in.intervals(ys,cutoffs)
  stratum.counts <- as.vector(tapply(in.s1,which.intervals$interval.number,sum))
  cutoffs.augmented <- c(-Inf,cutoffs,Inf)
  prob.allh.given.x <- matrix(NA,n,H)
  prob.h.given.x <- pnorm(cutoffs.augmented[which.intervals$interval.number+1],mean=mean,sd=sd) -
    pnorm(cutoffs.augmented[which.intervals$interval.number],mean=mean,sd=sd)
  for(h in c(1:H)){
    prob.allh.given.x[,h] <- pnorm(cutoffs.augmented[h+1],mean=mean,sd=sd) -
      pnorm(cutoffs.augmented[h],mean=mean,sd=sd)
  }
  pi.given.x <- as.vector(prob.allh.given.x %*% fh) + fSRS
  #prob1.given.x <- pnorm(-0.09268116,mean,sd)
  #prob2.given.x <- pnorm(2.96097625,mean,sd) - pnorm(-0.09268116,mean,sd)
  #prob3.given.x <- 1 - pnorm(2.96097625,mean,sd)
  #prob.h.given.x <- prob1.given.x * (ys<-0.09268116) +
  #  prob2.given.x * ((ys>=-0.09268116)&(ys<2.96097625)) +
  #  prob3.given.x * (ys>=2.96097625)
  #mean.pi.given.x <- 82/318*prob1.given.x + 137/1368*prob2.given.x + 81/314*prob3.given.x
  if(type=="SL"){
    ls <- sum( dnorm(ys,mean,sd,log=T)  - log(pi.given.x) )
  }
  if(type=="CL"){
    ls <- sum( dnorm(ys,mean,sd,log=T) - in.s1*log(prob.h.given.x) )
  }
  if(log) return(ls) else return(exp(ls))
}

