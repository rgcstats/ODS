#' Calculates the estimated likelihood for a regression model under a general sample design
#' assuming one continuous and one binary covariate with binary Y
#' and Y|x following a logistic regression model.
#'
#' This function calculates the likelihood.
#'
#' @details
#' Add some details later.
#'
#' @param beta vector of regression coefficients for Y|x1,x2
#' @param prob.x2 the parameter of x2's assumed bernoulli distribution
#' @param mu0.x1 the expected value of x1 conditional on x2=0
#' @param mu1.x1 the expected value of x1 conditional on x2=1
#' @param sd.x1 the SD of x1 conditional on x2
#' @param ys vector of sample values of the dependent variable
#' @param x1s sample values for covariate 1 (continuous)
#' @param x2s sample values for covariate 2 (binary)
#' @param pi1 The probability of selection when Y=1
#' @param pi0 The probability of selection when Y=0
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
#' L.2x.poisson.logistic(beta=c(1,0.5,-0.5),prob.x2=0.45,mu0.x1=0,mu1.x1=0,sd.x1=1,
#'        ys=1*(sample.example$y>1.5),x1s=sample.example$x1,x2s=sample.example$x2,
#'        pi1=0.5,pi0=0.1,log=TRUE)
#' @export
L.2x.poisson.logistic <- function(beta,prob.x2,mu0.x1,mu1.x1,sd.x1,
                                  ys,x1s,x2s,pi1,pi0,log=FALSE,
                                  verbose=FALSE){
  n <- length(ys)
  mu <- make.link("logit")$linkinv(as.vector(beta[1]+beta[2]*x1s+beta[3]*x2s))
  l.misspecified <- log(prob.x2)*sum(x2s) + log(1-prob.x2)*sum(1-x2s) +
    sum( dnorm( x=x1s , mean=mu0.x1*(1-x2s)+mu1.x1*x2s , sd=sd.x1 , log=TRUE ) ) +
    sum(log(mu[ys==1])) + sum(log(1-mu[ys==0]))
  # integrate pi.fn over x1, x2, y
  # Firstly get P[Y=1|x2=0], by integrating f(beta1+beta2*x1+beta2*0) over
  #     x1~N(mu0.x1,sd.x1) conditional on x2
  integrand0 <- function(x1){
    make.link("logit")$linkinv(beta[1]+beta[2]*x1) * dnorm(x1,mu0.x1,sd.x1)
  }
  mu0.y <- integrate(integrand0,-Inf,Inf)$value
  # Similarly get P[Y=1|x2=1], by integrating f(beta1+beta2*x1+beta2*1) over
  #     x1~N(mu1.x1,sd.x1) conditional on x2
  integrand1 <- function(x1){
    make.link("logit")$linkinv(beta[1]+beta[2]*x1+beta[3]) * dnorm(x1,mu1.x1,sd.x1)
  }
  mu1.y <- integrate(integrand1,-Inf,Inf)$value
  # Now we can get pi.bar = E[q(Y)]
  pi.bar <- ( (1-prob.x2)*mu0.y + prob.x2*mu1.y )*pi1  +
    ( (1-prob.x2)*(1-mu0.y) + prob.x2*(1-mu1.y) )*pi0
  # Now finish up
  l <- l.misspecified - n*log(pi.bar)
  if(verbose) cat("l.misspecified=",l.misspecified," ; l=",l," beta=(",paste(beta,sep=""),") ; sd=",sd,"\n")
  if(log) return(l)
  if(!log) return(exp(l))
}
