% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/L_2x.R
\name{L.2x}
\alias{L.2x}
\title{Calculates the estimated likelihood for a regression model under a general sample design
assuming one continuous and one binary covariate}
\usage{
L.2x(beta, sd, prob.x2, mu0.x1, mu1.x1, sd.x1, ys, x1s, x2s, N, p.s,
  pi.fn = NULL, specs = NULL, log = FALSE, pi.s, R = 100, all.errors.y,
  all.unif.x2, all.subsamp.unif, all.errors.x1, verbose = FALSE)
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

\item{N}{the population size}

\item{p.s}{A function defining the ODS sample design.
It must take arguments ys, yr, log and stuff, and return
the probability (or log of the probability) that a particular
sample was selected.}

\item{pi.fn}{Optional. The approximate probability that a unit is selected as a function of y. If provided,
a more efficient algorithm is used.}

\item{specs}{An object containing detailed specifications of the design.}

\item{log}{If FALSE, the function returns the likelihood, otherwise the log-likelihood}

\item{pi.s}{the probabilities of selection of the sample units}

\item{R}{the number of simulations used in the monte carlo approximation
of the likelihood}

\item{all.errors.y}{Optional. An R by (N-n)*10 real-valued matrix
of standard normal realisations.}

\item{all.unif.x2}{Optional. An R by (N-n)*10 real-valued matrix of U(0,1) realisations}

\item{all.subsamp.unif}{Optional. An R by (N-n)*10 real-valued matrix of U(0,1) realisations
of standard normal realisations.}

\item{all.errors.x1}{Optional. An R by (N-n)*10 real-valued matrix}

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
L.2x(beta=c(1,0.5,-0.5),sd=1,prob.x2=0.45,mu0.x1=0,mu1.x1=0,sd.x1=1,
     ys=sample.example$y,x1s=sample.example$x1,x2s=sample.example$x2,N=1000,
     p.s=p.s.ushaped,specs=c(n0=200,tuner=0.1,return.what=1),
     pi.fn=eg.pi.fn,
     pi.s=sample.example$pi,R=10,log=TRUE)
}
