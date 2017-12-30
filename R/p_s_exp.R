#' Calculates p(s) for a Poisson design with pi proportional to exp(c1+c2y) for given c1 and c2
#'
#' This function returns the probability p(s) of an entire sample being selected
#' for the Poisson sampling design where probabilities are proportional to exp(c1+c2y).
#' Probabilities of selection are truncated to be less than or equal to 1.
#'
#'  @details
#'  Note that expected sample size is not controlled.
#'
#' @param ys vector of the sample values of the dependent variable
#' @param yr vector of the non-sample values of the dependent variable
#' @param log If FALSE (the default), returns p(s). If TRUE, log(p(s)) is returned.
#' @param specs An object containing detailed specifications of the design. In this design
#'        specs should be a 2-vector with named elements c1 and c2
#' @export
p.s.exp <- function(ys,yr,log=FALSE,specs=NULL){
  yu <- c(ys,yr)
  N <- length(yu)
  n <- length(ys)
  pi.u <- pmin( exp(specs["c1"]+specs["c2"]*yu) , 1 )
  log.p.s <- sum(log(pi.u[1:n])) + sum(log(1-pi.u[(n+1):N])) #+ lchoose(N,n)
  if(!log) return( exp(log.p.s) )
  if(log) return( log.p.s )
}
