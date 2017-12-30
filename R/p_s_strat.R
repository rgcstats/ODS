#' Calculates p(s) for the design stratified by y
#'
#' This function returns the probability p(s) of an entire sample being selected
#' for stratified simple random sampling without replacement where strata are
#' defined by dividing y into intervals.
#'
#'  @details
#'  Add some details later.
#'
#' @param ys vector of the sample values of the dependent variable
#' @param yr vector of the non-sample values of the dependent variable
#' @param log If FALSE (the default), returns p(s). If TRUE, log(p(s)) is returned.
#' @param specs An object containing detailed specifications of the design. Should be vector of
#' cutoffs (for H strata there should be a vector of H-1 cutoffs)
#' @return The probability that s was selected (or the log thereof) multiplied by a constant.
#' @export
p.s.strat <- function(ys,yr,log=F,specs){
  n <- length(ys)
  yu <- c(ys,yr)
  cutoffs <- c(-Inf,specs,Inf)
  h.u <- cut(yu,breaks=cutoffs,include.lowest=TRUE)
  h.s <- cut(ys,breaks=cutoffs,include.lowest=TRUE)
  Nh <- as.vector(table(h.u))
  nh <- as.vector(table(h.s))
  log.p.s <- -sum(lchoose(Nh,nh))
  if(log) return(log.p.s)
  if(!log) return(exp(log.p.s))
}
