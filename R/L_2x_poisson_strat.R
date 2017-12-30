#' Calculates the estimated likelihood for a regression model under a general sample design
#' assuming one continuous and one binary covariate
#'
#' This function calculates the likelihood.
#'
#' @details
#' Add some details later.
#'
#' @param beta vector of regression coefficients
#' @param sd error standard deviation
#' @param prob.x2 the parameter of x2's assumed bernoulli distribution
#' @param mu0.x1 the expected value of x1 conditional on x2=0
#' @param mu1.x1 the expected value of x1 conditional on x2=1
#' @param sd.x1 the SD of x1 conditional on x2
#' @param ys vector of sample values of the dependent variable
#' @param x1s sample values for covariate 1 (continuous)
#' @param x2s sample values for covariate 2 (binary)
#' @param pi.h selection probability in each stratum (H-vector)
#' @param cutoffs vector of H-1 cutoffs on Y defining the strata
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @param verbose If TRUE, the function outputs information to the console.
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
L.2x.poisson.strat <- function(beta,sd,prob.x2,mu0.x1,mu1.x1,sd.x1,
                         ys,x1s,x2s,pi.h,cutoffs,log=FALSE,verbose=FALSE){
  n <- length(ys)
  l.misspecified <- log(prob.x2)*sum(x2s) + log(1-prob.x2)*sum(1-x2s) +
    sum( dnorm( x=x1s , mean=mu0.x1*(1-x2s)+mu1.x1*x2s , sd=sd.x1 , log=TRUE ) ) +
    sum( dnorm(x=ys,mean=as.vector(beta[1]+beta[2]*x1s+beta[3]*x2s),sd=sd,log=TRUE) )
  # integrate pi.fn over x1, x2, y (1000 simulations?)
  #   this should be easy because (y,x1) are bivariate normal given x2, and x2 is binary!
  mu0.y <- beta[1]+beta[2]*mu0.x1 # E[Y|x2=0]
  mu1.y <- beta[1]+beta[2]*mu1.x1+beta[3] # E[Y|x2=1]
  sd.condx2.y <- sqrt( beta[2]^2*sd.x1^2 + sd^2 ) # SD[Y|x2]
  # calculate probs of being less than or equal to each cutoff
  #   (appending +- inf)
  H <- length(cutoffs)+1
  cdf.cutoffs0 <- c( 0 , pnorm(cutoffs,mean=mu0.y,sd=sd.condx2.y) , 1 )
  cdf.cutoffs1 <- c( 0 , pnorm(cutoffs,mean=mu1.y,sd=sd.condx2.y) , 1 )
  prob.strata0 <- diff(cdf.cutoffs0)
  prob.strata1 <- diff(cdf.cutoffs1)
  prob.strata <- prob.x2*prob.strata1+(1-prob.x2)*prob.strata0
  pi.bar <- sum(prob.strata*pi.h)
  l <- l.misspecified - n*log(pi.bar)
  if(verbose) cat("l.misspecified=",l.misspecified," ; l=",l," beta=(",paste(beta,sep=""),") ; sd=",sd,"\n")
  if(log) return(l)
  if(!log) return(exp(l))
}
