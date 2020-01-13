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
#' @param pi.fn The approximate probability that a unit is selected as a function of y
#' @param specs An object containing detailed specifications of the design.
#' @param log If FALSE, the function returns the likelihood, otherwise the log-likelihood
#' @param pi.s the probabilities of selection of the sample units
#' @param R the number of simulations used in the monte carlo approximation
#'        of the likelihood
#' @param Rquad number of quadrature points
#' @param all.errors.y Optional. An R by (N-n)*10 real-valued matrix
#' of standard normal realisations.
#' @param all.unif.x2 Optional. An R by (N-n)*10 real-valued matrix of U(0,1) realisations
#' @param all.errors.x1 Optional. An R by (N-n)*10 real-valued matrix
#' @param all.subsamp.unif Optional. An R by (N-n)*10 real-valued matrix of U(0,1) realisations
#' of standard normal realisations.
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
#' L.2x.efficient(beta=c(1,0.5,-0.5),sd=1,prob.x2=0.45,mu0.x1=0,mu1.x1=0,sd.x1=1,
#'        ys=sample.example$y,x1s=sample.example$x1,x2s=sample.example$x2,N=1000,
#'        p.s=p.s.ushaped,specs=c(n0=200,tuner=0.1,return.what=1),
#'        pi.fn=eg.pi.fn,
#'        pi.s=sample.example$pi,R=10,log=TRUE)
#' @export
L.2x.efficient <- function(beta,sd,prob.x2,mu0.x1,mu1.x1,sd.x1,
                  ys,x1s,x2s,N,p.s,pi.fn,specs=NULL,log=FALSE,pi.s,R=100,Rquad=100,
              all.errors.y=NULL,all.unif.x2=NULL,
              all.subsamp.unif=NULL,all.errors.x1=NULL,verbose=FALSE){
  n <- length(ys)
  expand <- 4
  l.misspecified <- log(prob.x2)*sum(x2s) + log(1-prob.x2)*sum(1-x2s) +
    sum( dnorm( x=x1s , mean=mu0.x1*(1-x2s)+mu1.x1*x2s , sd=sd.x1 , log=TRUE ) ) +
    sum( dnorm(x=ys,mean=as.vector(beta[1]+beta[2]*x1s+beta[3]*x2s),sd=sd,log=TRUE) )
  # integrate pi.fn over x1, x2, y (1000 simulations?)
  #   this should be easy because (y,x1) are bivariate normal given x2, and x2 is binary!
  epsilon.vals <- qnorm(p=c(1:200)/201)
  mu0.y <- beta[1]+beta[2]*mu0.x1
  mu1.y <- beta[1]+beta[2]*mu1.x1+beta[3]
  sd.condx2.y <- sqrt( beta[2]^2*sd.x1^2 + sd^2 )
  # this is how pi.bar was calculated, new method is by quadrature with Rquad points
  #pi.bar <- (1-prob.x2)*mean(pmin(pi.fn(qnorm(p=c(1:40000)/40001,mean=mu0.y,sd=sd.condx2.y)),0.8)) +
  #  prob.x2*mean(pmin(pi.fn(qnorm(p=c(1:40000)/40001,mean=mu1.y,sd=sd.condx2.y)),0.8))
  rule <- fastGHQuad::gaussHermiteData(Rquad)
  ysim0 <- mu0.y+sd.condx2.y*sqrt(2)*rule$x
  pi.bar0 <- sum(pmin(pi.fn(ysim0),0.8)*rule$w)/sqrt(pi)
  ysim1 <- mu1.y+sd.condx2.y*sqrt(2)*rule$x
  pi.bar1 <- sum(pmin(pi.fn(ysim1),0.8)*rule$w)/sqrt(pi)
  pi.bar <- (1-prob.x2)*pi.bar0 + prob.x2*pi.bar1
  #
  # now numerically calculate difficult integral by monte carlo
  pi.s.guess <- pmin(pi.fn(ys),0.8)
  if(R==0) log.sampdesign.factor <- 0
  if(R>0){
    all.log.p.s <- all.log.p.s.guess <- rep(NA,R)
    if(is.null(all.subsamp.unif)) all.subsamp.unif <-
        matrix(runif((N-n)*10*R),nrow=R,ncol=(N-n)*expand)
    for(r in c(1:R)){
      # generate x2 first
      if(is.null(all.unif.x2)) x2r.sim <- 1*(runif((N-n)*expand)<=prob.x2)
      if(!is.null(all.unif.x2)) x2r.sim <- 1*(all.unif.x2[r,]<=prob.x2)
      # now generate x1|x2
      if(is.null(all.errors.x1)) x1r.sim <- mu0.x1*(1-x2r.sim) + mu1.x1*x2r.sim + sd.x1*rnorm((N-n)*expand)
      if(!is.null(all.errors.x1)) x1r.sim <- mu0.x1*(1-x2r.sim) + mu1.x1*x2r.sim + sd.x1*all.errors.x1[r,]
      if(is.null(all.errors.y)) errors.y <- rnorm(n=(N-n)*expand)
      if(!is.null(all.errors.y)) errors.y <- all.errors.y[r,]
      yr.sim <- beta[1] + beta[2]*x1r.sim + beta[3]*x2r.sim + errors.y*sd
      pi.r.guess <- pmin(pi.fn(yr.sim),0.8)
      yr.sim <- yr.sim[all.subsamp.unif<=(1-pi.r.guess)]
      yr.sim <- yr.sim[1:(N-n)]
      pi.r.guess <- pmin(pi.fn(yr.sim),0.8)
      all.log.p.s[r] <- p.s(ys=ys,yr=yr.sim,specs=specs,log=T)
      all.log.p.s.guess[r] <- sum(log(1-pi.r.guess)) + sum(log(pi.s.guess))
    }
    all.log.ratio <- all.log.p.s - all.log.p.s.guess
    if(all(is.infinite(all.log.ratio))) K <- 0 else K <- median(all.log.ratio[!is.infinite(all.log.ratio)])
    log.sampdesign.factor <- log( mean(exp(all.log.ratio-K)) ) + K
  }
  l <- l.misspecified + sum(log(pi.s.guess)) + (N-n)*log(1-pi.bar) + log.sampdesign.factor
  if(verbose) cat("l.misspecified=",l.misspecified," ; l=",l," beta=(",paste(beta,sep=""),") ; sd=",sd,"\n")
  if(log) return(l)
  if(!log) return(exp(l))
}
