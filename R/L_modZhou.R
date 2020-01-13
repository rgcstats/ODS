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
#'  If there are H cutoffs there are H+1 strata. NB ONLY H=3 IS ALLOWED!
#' @param R number of points used in quadrature calculation of integral
#' @param in.s0 A vector equal to 1 for units selected in the SRSWOR and to 0 for those
#'  in the stratified SRSWOR. Default is that there are no SRSWOR units.
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @return likelihood or log-likelihood conditional on stratum memberships
#' @examples
#' data(population_example)
#' L.modZhou(ys=sample.example$y,mean=0+0.5*sample.example$x,sd=1,log=TRUE,cutoffs=c(0,2))
#' @export
L.modZhou <- function(ys,mean,sd,cutoffs,R=100,in.s0,log=FALSE){
  n <- length(ys)
  H <- length(cutoffs)+1
  if(missing(in.s0)) in.s0 <- rep(0,n)
  in.s1 <- 1 - in.s0
  n0 <- sum(in.s0)
  n1 <- sum(in.s1)
  if(length(sd)==1) sd <- rep(sd,n)
  if(length(mean)==1) mean <- rep(mean,n)
  which.intervals <- locate.points.in.intervals(ys,cutoffs)
  stratum.counts <- as.vector(tapply(in.s1,which.intervals$interval.number,sum))
  mean.prob.strat.h.given.x <- matrix(NA,n,H)
  mean.prob.strat.h.given.x[,1] <- pnorm(cutoffs[1],mean=mean,sd=sd)
  if(H>2){
    for(h in c(2:(H-1))){
      mean.prob.strat.h.given.x[,h] <- pnorm(cutoffs[h],mean=mean,sd=sd) - pnorm(cutoffs[h-1],mean=mean,sd=sd)
    }
  }
  mean.prob.strat.h.given.x[,H] <- 1 - pnorm(cutoffs[H-1],mean=mean,sd=sd)
  objfn <- function(mean.prob.strat.h.except.last.stratum,return.p=FALSE){
    mean.prob.strat.h <- c(mean.prob.strat.h.except.last.stratum,1-sum(mean.prob.strat.h.except.last.stratum))
    p.i.inverse <- rep(n0,n)
    for(h in c(1:H)){
      p.i.inverse <- p.i.inverse + stratum.counts[h]*mean.prob.strat.h.given.x[,h]/mean.prob.strat.h[h]
    }
    p.i <- 1 / p.i.inverse
    discrepancy <- as.vector( t(mean.prob.strat.h.given.x) %*% p.i - mean.prob.strat.h )
    out <- sum(discrepancy^2)
    if(!return.p) return(out) else return(list(out=out,p.i=p.i,mean.prob.strat.h=mean.prob.strat.h))
  }
  #startval <- pmax(stratum.counts/n1,0.01)/sum(pmax(stratum.counts/n1,0.01))
  #lower <- pmax( 1e-5 , apply(mean.prob.strat.h.given.x,2,min)[1:(H-1)]*0.99 )
  #upper <- pmin( 1-1e-5 , apply(mean.prob.strat.h.given.x,2,max)[1:(H-1)] )
  #mean.prob.strat.h.except.last.stratum <- optim( par=startval[-H] , fn=objfn ,
  #                               method="L-BFGS-B",lower=lower,upper=upper )$par
  #
  # search 0.02*0.02 boxes progressively
  best.val <- Inf
  gridvals.lo <- c(0.01,seq(from=0.05,to=0.95,by=0.05))
  gridvals.hi <- c(0.05,seq(from=0.1,to=0.95,by=0.05),0.99)
  for(k1 in c(1:length(gridvals.lo))){
    for(k2 in c(1:length(gridvals.lo))){
      lower <- gridvals.lo[c(k1,k2)]
      upper <- gridvals.hi[c(k1,k2)]
      start <- c(lower+upper)/2
      if((start[1]+start[2])<1){
        o <- optim(par = start,
                   fn = objfn, method = "L-BFGS-B",
                   lower = lower , upper = upper )
        if(o$value<best.val){
          best.val <- o$value
          best.par <- o$par
          #cat("start1=",start1,"start2=",start2," objective=",o$value,"\n")
        }
      }
    }
  }
  mean.prob.strat.h.except.last.stratum <- best.par
  opt.results <- objfn(mean.prob.strat.h.except.last.stratum,return.p=TRUE)
  p.i <- opt.results$p.i
  mean.prob.strat.h <- opt.results$mean.prob.strat.h
  out <- sum( dnorm(x=ys,mean=mean,sd=sd,log=TRUE) ) +
    sum(log(p.i)) - sum(stratum.counts*log(mean.prob.strat.h))
  if(log) return(out) else return(exp(out))
}
