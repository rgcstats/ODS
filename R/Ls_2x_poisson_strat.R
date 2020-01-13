#' Calculates the estimated likelihood for a regression model under a general sample design
#' assuming one continuous and one binary covariate
#'
#' This function calculates the likelihood.
#'
#' @details
#' Add some details later.
#'
#' @param ys vector of sample values of the dependent variable
#' @param mean vector of means of dependent variable
#' @param sd vector of standard deviations of dependent variable
#' @param pi.h selection probability in each stratum (H-vector)
#' @param cutoffs vector of H-1 cutoffs on Y defining the strata
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @return The likelihood or log-likelihood.
#' @examples
#' data(population_example)
#' eg.ybar.pi <- sum(sample.example$y/sample.example$pi) / sum(1/sample.example$pi)
#' eg.s2.pi <- sum((sample.example$y-eg.ybar.pi)^2/sample.example$pi)/sum(1/sample.example$pi)
#' eg.pi.fn <- function(y){ out <- sqrt( 0.1 + 0.9*(y-eg.ybar.pi)^2/eg.s2.pi )
#'                          out <- out * mean(sample.example$pi) /
#'                                 mean(sqrt( 0.1 + 0.9*(sample.example$y-eg.ybar.pi)^2/eg.s2.pi ))
#'                          out  }
#' L.2x.poisson(beta=c(1,0.5,-0.5),sd=1,prob.x2=0.45,mu0.x1=0,mu1.x1=0,sd.x1=1,
#'        ys=sample.example$y,x1s=sample.example$x1,x2s=sample.example$x2,
#'        pi.fn=eg.pi.fn,pi.s=sample.example$pi,R=10,log=TRUE)
#' @export
Ls.2x.poisson.strat <- function(ys,mean,sd,pi.h,cutoffs,log=FALSE){
  n <- length(ys)
  H <- length(cutoffs)+1
  if(length(sd)==1) sd <- rep(sd,n)
  mean.pi.given.x <- rep(0,n)
  extended.cutoffs <- c(-Inf,cutoffs,Inf)
  for(h in c(1:H)){
    probs.h.given.x <- pnorm(q=extended.cutoffs[h+1],mean=mean,sd=sd)-
      pnorm(q=extended.cutoffs[h],mean=mean,sd=sd)
    mean.pi.given.x <- mean.pi.given.x + probs.h.given.x*pi.h[h]
  }
  ls <- sum( dnorm(ys,mean,sd,log=T) + log(pi.fn(ys)) - log(mean.pi.given.x) )
  if(log) return(ls) else return(exp(ls))
}
