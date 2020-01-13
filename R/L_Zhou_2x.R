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
#' @param cutoffs A vector of finite cutoffs defining the strata.
#'  If there are H cutoffs there are H+1 strata.
#' @param in.s0 A vector equal to 1 for units selected in the SRSWOR and to 0 for those
#'  in the stratified SRSWOR. Default is that there are no SRSWOR units.
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @return The likelihood or log-likelihood.
#' @examples
#' data(population_example)
#' eg.ybar.pi <- sum(sample.example$y/sample.example$pi) / sum(1/sample.example$pi)
#' eg.s2.pi <- sum((sample.example$y-eg.ybar.pi)^2/sample.example$pi)/sum(1/sample.example$pi)
#' L.Zhou.2x(beta=c(1,0.5,-0.5),sd=1,prob.x2=0.45,mu0.x1=0,mu1.x1=0,sd.x1=1,
#'        ys=sample.example$y,x1s=sample.example$x1,x2s=sample.example$x2,
#'        cutoffs=c(0,2),log=TRUE)
#' @export
L.Zhou.2x <- function(beta,sd,prob.x2,mu0.x1,mu1.x1,sd.x1,
                           ys,x1s,x2s,cutoffs,in.s0,log=FALSE){
  n <- length(ys)
  H <- length(cutoffs)+1
  #if(length(cutoffs)!=2) stop("This function requires exactly 3 strata.")
  if(missing(in.s0)) in.s0 <- rep(0,n)
  in.s1 <- 1 - in.s0
  which.intervals <- locate.points.in.intervals(ys,cutoffs)
  stratum.counts <- as.vector(tapply(in.s1,which.intervals$interval.number,sum))
  l.misspecified <- log(prob.x2)*sum(x2s) + log(1-prob.x2)*sum(1-x2s) +
    sum( dnorm( x=x1s , mean=mu0.x1*(1-x2s)+mu1.x1*x2s , sd=sd.x1 , log=TRUE ) ) +
    sum( dnorm(x=ys,mean=as.vector(beta[1]+beta[2]*x1s+beta[3]*x2s),sd=sd,log=TRUE) )
  mu0.y <- beta[1]+beta[2]*mu0.x1 # E[Y|x2=0]
  mu1.y <- beta[1]+beta[2]*mu1.x1+beta[3] # E[Y|x2=1]
  sd.condx2.y <- sqrt( beta[2]^2*sd.x1^2 + sd^2 ) # SD[Y|x2]
  prob.le.cutoffs <- pnorm(cutoffs,mean=mu0.y,sd=sd.condx2.y)*(1-prob.x2) +
    pnorm(cutoffs,mean=mu1.y,sd=sd.condx2.y)*prob.x2
  prob.le.cutoffs.augmented <- c(0,prob.le.cutoffs,1)
  prob.strat.h <- prob.le.cutoffs.augmented[-1] - prob.le.cutoffs.augmented[-length(prob.le.cutoffs.augmented)]
  #prob.strat1 <- pnorm(cutoffs[1],mean=mu0.y,sd=sd.condx2.y)*(1-prob.x2) +
  #  pnorm(cutoffs[1],mean=mu1.y,sd=sd.condx2.y)*prob.x2
  #prob.strat3 <- pnorm(cutoffs[2],mean=mu0.y,sd=sd.condx2.y,lower.tail=FALSE)*(1-prob.x2) +
  #  pnorm(cutoffs[2],mean=mu1.y,sd=sd.condx2.y,lower.tail=FALSE)*prob.x2
  #prob.strat2 <- 1 - prob.strat1 - prob.strat3
  out <- l.misspecified - sum(stratum.counts*log(prob.strat.h))
    #stratum.counts[1]*log(prob.strat1) - stratum.counts[2]*log(prob.strat2) -
    #stratum.counts[3]*log(prob.strat3)
  if(log) return(out) else return(exp(out))
}
