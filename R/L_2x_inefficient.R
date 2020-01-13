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
#' @param N the population size
#' @param p.s A function defining the ODS sample design.
#' It must take arguments ys, yr, log and stuff, and return
#' the probability (or log of the probability) that a particular
#' sample was selected.
#' @param specs An object containing detailed specifications of the design.
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @param pi.s the probabilities of selection of the sample units
#' @param pi.fn not used
#' @param R the number of simulations used in the monte carlo approximation
#'        of the likelihood
#' @param all.errors.y Optional. An R by N-n real-valued matrix
#' of standard normal realisations.
#' @param all.unif.x2 Optional. An R by N-n real-valued matrix of U(0,1) realisations
#' @param all.subsamp.unif not used.
#' @param all.errors.x1 Optional. An R by N-n real-valued matrix
#' of standard normal realisations.
#' @param verbose If TRUE, the function outputs information to the console.
#' @return The likelihood or log-likelihood.
#' @examples
#' data(population_example)
#' L.2x.inefficient(beta=c(1,0.5,-0.5),sd=1,prob.x2=0.45,mu0.x1=0,mu1.x1=0,sd.x1=1,
#'        ys=sample.example$y,x1s=sample.example$x1,x2s=sample.example$x2,N=1000,
#'        p.s=p.s.ushaped,pi.fn=NULL,specs=c(n0=200,tuner=0.1,return.what=1),
#'        pi.s=sample.example$pi,R=10,log=TRUE)
#' @export
L.2x.inefficient <- function(beta,sd,prob.x2,mu0.x1,mu1.x1,sd.x1,
                             ys,x1s,x2s,N,p.s,pi.fn=NULL,specs=NULL,log=FALSE,pi.s,R=100,
                             all.errors.y=NULL,all.unif.x2=NULL,
                             all.subsamp.unif=NULL,
                             all.errors.x1=NULL,verbose=FALSE){
  n <- length(ys)
  l.misspecified <- log(prob.x2)*sum(x2s) + log(1-prob.x2)*sum(1-x2s) +
    sum( dnorm( x=x1s , mean=mu0.x1*(1-x2s)+mu1.x1*x2s , sd=sd.x1 , log=TRUE ) ) +
    sum( dnorm(x=ys,mean=as.vector(beta[1]+beta[2]*x1s+beta[3]*x2s),sd=sd,log=TRUE) )
  all.log.p.s <- rep(NA,R)
  for(r in c(1:R)){
    # generate x2 first
    if(is.null(all.unif.x2)) x2r.sim <- 1*(runif(N-n)<=prob.x2)
    if(!is.null(all.unif.x2)) x2r.sim <- 1*(all.unif.x2[r,]<=prob.x2)
    # now generate x1|x2
    if(is.null(all.errors.x1)) x1r.sim <- mu0.x1*(1-x2r.sim) + mu1.x1*x2r.sim + sd.x1*rnorm(N-n)
    if(!is.null(all.errors.x1)) x1r.sim <- mu0.x1*(1-x2r.sim) + mu1.x1*x2r.sim + sd.x1*all.errors.x1[r,]
    if(is.null(all.errors.y)) errors.y <- rnorm(n=N-n)
    if(!is.null(all.errors.y)) errors.y <- all.errors.y[r,]
    yr.sim <- beta[1] + beta[2]*x1r.sim + beta[3]*x2r.sim + errors.y*sd
    all.log.p.s[r] <- p.s(ys=ys,yr=yr.sim,specs=specs,log=T)
  }
  K <- median(all.log.p.s)
  log.sampdesign.factor <- log( mean(exp(all.log.p.s-K)) ) + K
  l <- l.misspecified + log.sampdesign.factor
  if(verbose) cat("l.misspecified=",l.misspecified," ; l=",l," beta=(",paste(beta,sep=""),") ; sd=",sd,"\n")
  if(log) return(l)
  if(!log) return(exp(l))
}
