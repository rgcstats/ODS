#' Calculates p(s) for a Poisson design with pi proportional to exp(c2*y/SY) for given c2
#'
#' This function returns the probability p(s) of an entire sample being selected
#' for the Poisson sampling design where probabilities are proportional to exp(c2*y/SY).
#' SY is the population SD of Y, and constant of proportionality calculated such that
#' expected sample size is n0.
#' Probabilities of selection are truncated to be less than or equal to 1.
#'
#'  @details
#'  Add some details later.
#'
#' @param ys vector of the sample values of the dependent variable
#' @param yr vector of the non-sample values of the dependent variable
#' @param log If FALSE (the default), returns p(s). If TRUE, log(p(s)) is returned.
#' @param specs An object containing detailed specifications of the design. In this design
#'        specs should be a 3-vector with named elements c2, n0 and return.what
#' @return If specs$return.what=1, returns the probability that s was selected (or the log thereof).
#'         If specs$return.what=2, returns the probabilities of selection for each population unit
#'         If specs$return.what=3, returns a function pi.fn(y) which gives the unit probabilities of selection as a function of y
#' @export
p.s.exp.adaptive <- function(ys,yr,log=FALSE,specs=NULL){
  yu <- c(ys,yr)
  N <- length(yu)
  n <- length(ys)
  SY <- sqrt(var(yu))
  zu <- (yu-mean(yu))/sqrt(var(yu))
  c2 <- specs["c"]
  n0 <- specs["n0"]
  return.what <- specs["return.what"]
  totsamp.discrepancy <- function(lambda,size,n0,ub,lb){
    ( sum(pmax(pmin(lambda*n0*size/sum(size),ub),lb)) - n0 )^2
  }
  size <- exp(c2*zu)
  lambda <- optimise(f=totsamp.discrepancy,interval=c(1,10),size=size,n0=n0,ub=0,ub=1)$minimum
  pi.u <- pmin(lambda*n0*size/sum(size),1)
  if(return.what==2) return(pi.u)
  if(return.what==3) return(
    function(y){
      size.y <- exp(c2*(y-mean(yu))/sqrt(var(yu)))
      pmin(lambda*n0/sum(size)*size.y,1)
    }  )
  if(return.what==1){
    log.p.s <- sum(log(pi.u[1:n])) + sum(log(1-pi.u[(n+1):N])) #+ lchoose(N,n)
    if(!log) return( exp(log.p.s) )
    if(log) return( log.p.s )
  }
}
