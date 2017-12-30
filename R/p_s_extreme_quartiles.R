#' Calculates p(s) for the design selecting from quartiles 1 and 4
#'
#' This function returns the probability p(s) of an entire sample being selected
#' for the design with n/2 selected from first quartile and n/2 selected from fourth quartile
#' of the population values of the dependent variable
#'
#'  @details
#'  Insert some details here.
#'
#' @param ys vector of the sample values of the dependent variable
#' @param yr vector of the non-sample values of the dependent variable
#' @param log If FALSE (the default), returns p(s). If TRUE, log(p(s)) is returned.
#' @param specs An object containing detailed specifications of the design (not used for this design)
#' @return The probability that s was selected (or the log thereof) multiplied by N/4-choose-n/4.
#' (The multiplication by N/4-choose-n/4 is to avoid very small values of p(s).)
#' @export
p.s.extreme.quartiles <- function(ys,yr,log=F,specs=NULL){
  cuts <- quantile(c(ys,yr),probs=c(0.25,0.75))
  n <- length(ys)
  n1 <- sum(ys<cuts[1])
  n2 <- sum(ys>cuts[2])
  if((n1!=(n/2))|(n2!=(n/2))){
    if(!log) return(0)
    if(log) return(-Inf)
  } else{
    N <- n + length(yr)
    if(log) return(-lchoose(N/4,n/2)-lchoose(N/4,n/2))
    if(!log) return(1/choose(N/4,n/2)/choose(N/4,n/2))
  }
}
