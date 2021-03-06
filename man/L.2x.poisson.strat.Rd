% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L_2x_poisson_strat.R
\name{L.2x.poisson.strat}
\alias{L.2x.poisson.strat}
\title{Calculates the estimated likelihood for a regression model under a general sample design
assuming one continuous and one binary covariate}
\usage{
L.2x.poisson.strat(beta, sd, prob.x2, mu0.x1, mu1.x1, sd.x1, ys, x1s, x2s, pi.h,
  cutoffs, log = FALSE, verbose = FALSE)
}
\arguments{
\item{beta}{vector of regression coefficients}

\item{sd}{error standard deviation}

\item{prob.x2}{the parameter of x2's assumed bernoulli distribution}

\item{mu0.x1}{the expected value of x1 conditional on x2=0}

\item{mu1.x1}{the expected value of x1 conditional on x2=1}

\item{sd.x1}{the SD of x1 conditional on x2}

\item{ys}{vector of sample values of the dependent variable}

\item{x1s}{sample values for covariate 1 (continuous)}

\item{x2s}{sample values for covariate 2 (binary)}

\item{pi.h}{selection probability in each stratum (H-vector)}

\item{cutoffs}{vector of H-1 cutoffs on Y defining the strata}

\item{log}{If FALSE, the function returns the likelihood, otherwise the log-likelihood}

\item{verbose}{If TRUE, the function outputs information to the console.}
}
\value{
The likelihood or log-likelihood.
}
\description{
This function calculates the likelihood.
}
\details{
Add some details later.
}
\examples{
data(population_example)
eg.ybar.pi <- sum(sample.example$y/sample.example$pi) / sum(1/sample.example$pi)
eg.s2.pi <- sum((sample.example$y-eg.ybar.pi)^2/sample.example$pi)/sum(1/sample.example$pi)
eg.pi.fn <- function(y){ out <- sqrt( 0.1 + 0.9*(y-eg.ybar.pi)^2/eg.s2.pi )
                         out <- out * mean(sample.example$pi) /
                                mean(sqrt( 0.1 + 0.9*(sample.example$y-eg.ybar.pi)^2/eg.s2.pi ))
                         out  }
L.2x.poisson(beta=c(1,0.5,-0.5),sd=1,prob.x2=0.45,mu0.x1=0,mu1.x1=0,sd.x1=1,
       ys=sample.example$y,x1s=sample.example$x1,x2s=sample.example$x2,
       pi.fn=eg.pi.fn,pi.s=sample.example$pi,R=10,log=TRUE)
}
