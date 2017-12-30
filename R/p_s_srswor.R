#' Calculates p(s) for simple random sampling without replacement (SRSWOR)
#'
#' This function returns the probability p(s) of an entire sample being selected
#' under SRSWOR
#'
#'  @details
#'  Add some details later.
#'
#' @param ys vector of the sample values of the dependent variable
#' @param yr vector of the non-sample values of the dependent variable
#' @param log If FALSE (the default), returns p(s). If TRUE, log(p(s)) is returned.
#' @param specs An object containing detailed specifications of the design. Should contain n/N only.
#' @return The probability that s was selected (or the log thereof) multiplied by a constant.
#' Because p(s) is equal either to a constant or to zero, the function just returns a value of p(s)=1 or 0.
#' @export
p.s.srswor <- function(ys,yr,log=F,specs){
  if(!log) return(specs)
  if(log) return(log(specs))
  specs
}
