#' Calculates p(s) for the design stratified by quartile
#'
#' This function returns the probability p(s) of an entire sample being selected
#' for stratified simple random sampling without replacement where strata are
#' defined by the quartiles of the population values of the dependent variable
#'
#'  @details
#'  Add some details later.
#'
#' @param ys vector of the sample values of the dependent variable
#' @param yr vector of the non-sample values of the dependent variable
#' @param log If FALSE (the default), returns p(s). If TRUE, log(p(s)) is returned.
#' @param specs An object containing detailed specifications of the design. Should be a 4-vector
#'              containing the sample sizes in the quartiles 1-4.
#' @return The probability that s was selected (or the log thereof) multiplied by a constant.
#' @export
p.s.stratq <- function(ys,yr,log=F,specs){
  cuts <- quantile(c(ys,yr),probs=c(0.25,0.5,0.75))
  n <- length(ys)
  observed.nh <- c( sum(ys<cuts[1]) , sum((ys>=cuts[1])&(ys<cuts[2])) ,
                    sum((ys>=cuts[2])&(ys<cuts[3])) , sum(ys>=cuts[3]) )
  Nh <-  c( sum(yu<cuts[1]) , sum((yu>=cuts[1])&(yu<cuts[2])) ,
            sum((yu>=cuts[2])&(yu<cuts[3])) , sum(yu>=cuts[4]) )
  if(any(observed.nh!=specs)){
    if(!log) return(0)
    if(log) return(-Inf)
  } else{
    if(!log) return(1)
    if(log) return(0)
  }
}
