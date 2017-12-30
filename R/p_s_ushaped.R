#' Calculates p(s) for a u-shaped design
#'
#' This function returns the probability p(s) of an entire sample being selected
#' for the Poisson sampling design where probabilities are a u-shaped function of
#' y scaled such that the expected sample size is as specified.
#'
#'  @details
#' probability is proportional to sqrt(tuner+(1-tuner)*(y-meany.u)^2/vary.u)
#' truncated to lie be between 0.1*n0/N and 1-0.1*(N-n0)/N,
#'     or equivalently f0/10 and (f0+9)/10, where f0=n0/N,
#'
#' @param ys vector of the sample values of the dependent variable
#' @param yr vector of the non-sample values of the dependent variable
#' @param log If FALSE (the default), returns p(s). If TRUE, log(p(s)) is returned.
#' @param specs An object containing detailed specifications of the design. In this design
#'        specs should be a 3-vector with named elements n0, tuner and return.what.
#' @return If specs$return.what=1, returns the probability that s was selected (or the log thereof).
#'         If specs$return.what=2, returns the probabilities of selection for each population unit
#'         If specs$return.what=3, returns a function pi.fn(y) which gives the unit probabilities of selection as a function of y
#' @export
p.s.ushaped <- function(ys,yr,log=FALSE,specs=NULL){
  yu <- c(ys,yr)
  N <- length(yu)
  n0 <- specs["n0"]
  tuner <- specs["tuner"]
  return.what <- specs["return.what"]
  size <- sqrt(tuner+(1-tuner)*(yu-mean(yu))^2/var(yu))
  f0 <- n0/N
  lb <- f0/10
  ub <- (f0+9)/10
  totsamp.discrepancy <- function(lambda,size,n0){
    ( sum(pmax(pmin(lambda*n0*size/sum(size),ub),lb)) - n0 )^2
  }
  lambda <- optimise(f=totsamp.discrepancy,interval=c(1,10),size=size,n0=n0)$minimum
  pi.u <- pmax(pmin(lambda*n0*size/sum(size),ub),lb)
  if(return.what==2) return(pi.u)
  if(return.what==3) return(
    function(y){
      size <- sqrt(tuner+(1-tuner)*(y-mean(yu))^2/var(yu))
      pmax(pmin(lambda*n0*size/sum(size),ub),lb)
    }  )
  if(return.what==1){
    n <- length(ys)
    log.p.s <- sum(log(pi.u[1:n])) + sum(log(1-pi.u[(n+1):N])) #+ lchoose(N,n)
    if(!log) return( exp(log.p.s) )
    if(log) return( log.p.s )
  }
}
