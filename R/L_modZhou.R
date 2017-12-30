#' Calculates a modified Zhou conditional likelihood
#'
#' This function calculates a conditional likelihood
#' which is a modification of the one in Zhou et al. (2002).
#'
#' @details
#' The dependent variable is assumed to be normally distributed
#'
#' @param ys vector of sample values of the dependent variable
#' @param mean vector of means of dependent variable
#' @param sd vector of standard deviations of dependent variable
#' @param cutoffs A vector of finite cutoffs defining the strata.
#'  If there are H cutoffs there are H+1 strata.
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @return likelihood or log-likelihood conditional on stratum memberships
#' @examples
#' data(population_example)
#' L.modZhou(ys=sample.example$y,mean=5+0.5*sample.example$x,sd=1,log=TRUE,cutoffs=c(1,2,3))
#' @export
L.modZhou <- function(ys,mean,sd,cutoffs,log=FALSE){
  n <- length(ys)
  if(length(sd)==1) sd <- rep(sd,n)
  if(length(mean)==1) mean <- rep(mean,n)
  which.intervals <- locate.points.in.intervals(ys,cutoffs)
  # out <- sum( dtruncnorm(x=ys,a=lower.cutoffs,b=upper.cutoffs,mean=mean,sd=sd) )
  out <- sum( dnorm(x=ys,mean=mean,sd=sd,log=TRUE) -
                log(pnorm(which.intervals$upper.cutoffs)-pnorm(which.intervals$lower.cutoffs)) )
  if(log) return(out) else return(exp(out))
}
