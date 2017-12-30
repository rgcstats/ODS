#' Calculates p(s) for the biggest and smallest design
#'
#' This function returns the probability p(s) of an entire sample being selected
#'  for the design where the largest n/2 and smallest n/2 values of the
#'  dependent variable are selected from the population
#'
#'  @details
#'  Insert some details here.
#'
#' @param ys vector of the sample values of the dependent variable
#' @param yr vector of the non-sample values of the dependent variable
#' @param log If FALSE (the default), returns p(s). If TRUE, log(p(s)) is returned.
#' @param specs An object containing detailed specifications of the design (not used for this design)
#' @return The probability that s was selected (or the log thereof). Either 1 or 0 for this design.
#' @export
p.s.extremes <- function(ys,yr,log=FALSE,specs=NULL){
  if(any(ys==min(yr))|any(ys==max(yr))) stop("Error: duplicates affecting sampling")
  n <- length(ys)
  n1 <- sum(ys<min(yr))
  n2 <- sum(ys>max(yr))
  p.s <- 1*(n1==(n/2))*(n2==(n/2))
  if(log) return(log(p.s)) else return(p.s)
}
